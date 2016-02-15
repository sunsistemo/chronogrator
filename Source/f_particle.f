*     *==force.spg  processed by SPAG 4.52O  at 15:46 on 28 Mar 1996
      SUBROUTINE f_particle(i, Fx, Fy, Fz, En, Vir, nlist, list, swiver)
c
c     Calculate the force acting on particle i
c
c     Fx  (output) array: x component of the force acting on the particles
c     Fy  (output) array: y component of the force acting on the particles
c     Fz  (output) array: z component of the force acting on the particles
c     En  (output)      : total energy
c     Vir (output)      : total virial
c
c
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'system.inc'
      INCLUDE 'potential.inc'
      INCLUDE 'verlet.inc'
      INCLUDE 'samp.inc'

      DOUBLE PRECISION xi, yi, zi, En, dx, dy, dz, r2, Vir, virij, enij,
     &     fr, Fx, Fy, Fz, r2i, r6i
      INTEGER i, j, jj, swiver
      INTEGER nlist(npmax), list(npmax, npmax)
      DIMENSION Fx(*), Fy(*), Fz(*)

c     --- make forces zero
      Fx(i) = 0
      Fy(i) = 0
      Fz(i) = 0
c     --- Check whether to make new Verlet list
      IF (swiver.EQ.1) THEN
         IF (abs(X(i) - XV(i)).GT.(rdv)) THEN
            CALL VLIST(nlist, list, rv, swiver)
         END IF
      ELSE
         IF (abs(X(i) - XV2(i)).GT.(rdv2)) THEN
            CALL VLIST(nlist, list, rv2, swiver)
         END IF
      END IF
      xi = X(i)
      yi = Y(i)
      zi = Z(i)
c     --- For particle i calculate interaction with other particles
c     --- in its Verlet-list
      DO jj = 1, NLIST(i)
         j = list(i, jj)
         dx = xi - X(j)
         dy = yi - Y(j)
         dz = zi - Z(j)
c     ---periodic boundary conditions
         dx = dx - BOX*ANINT(dx/BOX)
         dy = dy - BOX*ANINT(dy/BOX)
         dz = dz - BOX*ANINT(dz/BOX)
         r2 = dx*dx + dy*dy + dz*dz
         IF (swiver.EQ.1 .AND. r2.LE.RC2) THEN
            r2i = 1/r2
            r6i = r2i*r2i*r2i
            enij = 4*(r6i*r6i-r6i) - ECUT
            virij = 48*(r6i*r6i-0.5D0*r6i)
            En = En + enij
            Vir = Vir + virij
            fr = virij*r2i
            Fx(i) = Fx(i) + fr*dx
            Fy(i) = Fy(i) + fr*dy
            Fz(i) = Fz(i) + fr*dz
            Fx(j) = Fx(j) - fr*dx
            Fy(j) = Fy(j) - fr*dy
            Fz(j) = Fz(j) - fr*dz
         END IF
      END DO
      RETURN
      END
