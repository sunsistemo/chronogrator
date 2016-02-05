c verlet.inc
      double precision xv,yv,zv,rv,rdv
      common/verlet1/xv(npmax),yv(npmax),zv(npmax),rv,rdv
c
c     xv(i),yv(i),zv(i)    : Verlet position particle i
c     rv                   : Verlet list radius
