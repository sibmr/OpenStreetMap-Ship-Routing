longsteps=18
latsteps = 9

longmin = -180
longmax = 180
latmin = -90
latmax = 90

longstep = (longmax-longmin)/longsteps
latstep = (latmax-latmin)/latsteps

print("[-180,-90]")
step_back = False
for i in range(longsteps):
    for j in range(latsteps):
        print(",")
        print(f"[{longmin+longstep*i},{latmin+latstep*(j+1)}]")
        print(",")
        print(f"[{longmin+longstep*(i+1)},{latmin+latstep*(j+1)}]")
        print(",")
        print(f"[{longmin+longstep*(i+1)},{latmin+latstep*(j)}]")
        print(",")
        print(f"[{longmin+longstep*(i)},{latmin+latstep*j}]")
        print(",")
        print(f"[{longmin+longstep*i},{latmin+latstep*(j+1)}]")
    print(",")
    print(f"[{longmin+longstep*i},{latmin+latstep}]")