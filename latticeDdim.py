#!/usr/bin/env python
"""
A quick exploration of higher dimensional hypercubic lattices

"""

def main():
     
    L = 4
    D = 6

    N = L**D
    print 'number of sites',N

    for i in range(N): #loop over all sites
        print i,
        for j in range(D): #loop over dimension
            if ( (i+L**j)%(L**(j+1)) >= 0) and ( (i+L**j)%(L**(j+1)) < L**j):
                print i+L**j - L**(j+1),
            else:
                print i+L**j,
        print 




if __name__ == '__main__':
    main()
