##
## patchy DLA model in Julia
## mjsottile@me.com
##


#
# create an n-by-n world with some initial value for all cells
# and an initial value for the central seed cell.
#
function setup(n, initval, seedval)
    x = fill(initval, (n,n))

    center = Int(ceil(n/2))
    x[center,center] = seedval

    x
end

#
# generate a random angle
#
rand_angle() = rand() * 2 * pi

#
# generate a point on a circle of radius r
#
function circpoint(r)
    theta = rand_angle()
    r * cos(theta) , r * sin(theta)
end

# return a boolean indicating if the given angle is inside the patch
# pointing in the orientation direction
function test_inside_patch(phi, orientation, patchsize)
    angle_lo = mod(orientation - patchsize/2, 2pi)
    angle_hi = mod(orientation + patchsize/2, 2pi)
    if angle_lo > angle_hi
        (0.0 <= phi <= angle_hi) || (angle_lo <= phi <= 2pi )
    else
        angle_lo <= phi <= angle_hi
    end
end

#
# collision test.  given the state of the world x, the point
# to test for collision, the angle of orientation for the
# patchy particle, and the patch size for the model
#
function collide(x, pt, angle, patchsize)
    # create a closure over the patchsize since that is the same for every
    # call.
    test_inside(phi,orientation) = 
        test_inside_patch(phi, orientation, patchsize)

    # look left: angle pi in patch?
    if pt[1] > 1 && x[pt[1]-1,pt[2]] < Inf
        test_inside(pi, angle) && test_inside(0.0, x[pt[1]-1,pt[2]])

    # look right: angle 0 in patch?
    elseif pt[1] < size(x)[1] && x[pt[1]+1,pt[2]] < Inf
        test_inside(0, angle) && test_inside(pi, x[pt[1]+1,pt[2]])

    # look up: angle pi/2 in patch?
    elseif pt[2] > 1 && x[pt[1],pt[2]-1] < Inf
        test_inside(pi/2, angle) && test_inside(3pi/2, x[pt[1],pt[2]-1])

    # look down: angle 3pi/2 in patch?
    elseif pt[2] < size(x)[2] && x[pt[1],pt[2]+1] < Inf
        test_inside(3pi/2, angle) && test_inside(pi/2, x[pt[1],pt[2]+1])

    else

        false
    end

end

#
# let the particle take a walk from its starting point to
# either collision, wandering off the edge, or exceeding
# the maximum number of steps allowed in the walk
#
function walk_particle(x, radius, angle, patchsize, maxiters)
    sz = size(x)

    # generate the point, shifting it to be inside the
    # index set for the array
    px,py = circpoint(radius)
    px,py = Int(ceil(px+sz[1]/2)), Int(ceil(py+sz[2]/2))

    iter = 0
    while iter < maxiters
        if collide(x, (px,py), angle, patchsize)
            x[px,py] = angle
            return (x,true)
        end

        dir = rand(1:4)

        xoffs = [1,-1,0,0]
        yoffs = [0,0,-1,1]

        px,py = px+xoffs[dir], py+yoffs[dir]

        if px < 1 || px > sz[1] || py < 1 || py > sz[2]
            return (x, false)
        end

        # oops - occupied cell, so exclusion principle kicks in.
        # this occurs if the cell is occupied but has a patch that
        # isn't facing the particle (or vice versa), so no stick.
        # back up, wasting a step.
        if x[px,py] < Inf
            px, py = px-xoffs[dir],py-yoffs[dir]
        end

        iter += 1
    end

    return (x, false)
end

# calculate for every cell how many potential attachment points there are
# given its neighbors
function compute_availables(x, patchsize)
    availables = zeros(size(x))

    for i=1:size(x)[1]
        for j=1:size(x)[2]
            # only empty cells are available
            if x[i,j] == Inf
                if i > 1 && test_inside_patch(0, x[i-1,j], patchsize)
                    availables[i,j] += 1
                end
                if i < size(x)[1] && test_inside_patch(pi, x[i+1,j], patchsize)
                    availables[i,j] += 1
                end
                if j > 1 && test_inside_patch(3pi/2, x[i,j-1], patchsize)
                    availables[i,j] += 1
                end
                if j < size(x)[2] && test_inside_patch(pi/2, x[i,j+1], patchsize)
                    availables[i,j] += 1
                end
            end
        end
    end

    availables
end

#
# main program
#
function main()
    n = 301

    seed_angle = rand_angle()

    x = setup(n, Inf, seed_angle)

    maxpart = 3000
    npart = 1
    patchsize = 1.99 * pi
    maxsteps = 100000

    maxattempts = maxpart * 3

    attempt = 1

    while npart < maxpart && attempt < maxattempts
        (x, hit) = walk_particle(x, Int(ceil(n*0.4)), rand_angle(), patchsize, maxsteps)
        if hit
            npart += 1
        end
        attempt += 1
    end

    # for i=1:size(x)[1]
    #     for j=1:size(x)[2]
    #         if x[i,j] < Inf
    #             println("$(x[i,j]) $i $j")
    #         end
    #     end
    # end
    available = compute_availables(x, patchsize)

    for i=1:size(x)[1]
        for j=1:size(x)[2]
            print("$(available[i,j]) ")
        end
        println()
    end

end

main()

