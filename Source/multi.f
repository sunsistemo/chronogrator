      SUBROUTINE MULTI(Fx, Fy, Fz, Fx2, Fy2, Fz2, nlist, list,
     &                 nlist2, list2, Enkin, Enpot, Vir, Delt, n)

      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'veloc.inc'
      INCLUDE 'verlet.inc'
      DOUBLE PRECISION Fx(*), Fy(*), Fz(*)
      DOUBLE PRECISION Fx2(*),Fy2(*), Fz2(*), En2, Vir2
      DOUBLE PRECISION Delt, dt, dt2, v2, Enkin, Enpot, Vir
      INTEGER nlist(npmax), list(npmax, npmax)
      INTEGER nlist2(npmax), list2(npmax, npmax)

      INTEGER i, j, n

      dt = Delt * n    ! big time step
      dt2 = Delt       ! small time step
      v2 = 0.D0
      Enkin = 0.D0

      DO i = 1, NPART
         VX(i) = VX(i) + 0.5 * dt * (Fx(i) - Fx2(i))
         VY(i) = VY(i) + 0.5 * dt * (Fy(i) - Fy2(i))
         VZ(i) = VZ(i) + 0.5 * dt * (Fz(i) - Fz2(i))
      END DO

      DO j = 1, n

         DO i = 1, NPART
            VX(i) = VX(i) + 0.5 * dt2 * Fx2(i)
            VY(i) = VY(i) + 0.5 * dt2 * Fy2(i)
            VZ(i) = VZ(i) + 0.5 * dt2 * Fz2(i)

            X(i) = X(i) + dt2 * 2 * VX(i)
            Y(i) = Y(i) + dt2 * 2 * VY(i)
            Z(i) = Z(i) + dt2 * 2 * VZ(i)
         END DO

         CALL FORCE(Fx2, Fy2, Fz2, En2, Vir2, nlist2, list2, .true., 2)

         DO i = 1, NPART
            VX(i) = VX(i) + 0.5 * dt2 * Fx2(i)
            VY(i) = VY(i) + 0.5 * dt2 * Fy2(i)
            VZ(i) = VZ(i) + 0.5 * dt2 * Fz2(i)
         END DO
      END DO

      CALL FORCE(Fx, Fy, Fz, Enpot, Vir, nlist, list, .true., 1)

      DO i = 1, NPART
         VX(i) = VX(i) + 0.5 * dt * (Fx(i) - Fx2(i))
         VY(i) = VY(i) + 0.5 * dt * (Fy(i) - Fy2(i))
         VZ(i) = VZ(i) + 0.5 * dt * (Fz(i) - Fz2(i))

         v2 = v2 + (VX(i)**2 + VY(i)**2 + VZ(i)**2)
      END DO
      Enkin = v2/2
      RETURN
      END
