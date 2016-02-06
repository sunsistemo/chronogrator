      SUBROUTINE MULTI(Fx, Fy, Fz, Fx2, Fy2, Fz2, nlist, list,
      &                nlist2, list2, Enkin, Delt, n)

      IMPLICIT NONE
      INCLUDE 'conf.inc'
      INCLUDE 'veloc.inc'
      INCLUDE 'parameter.inc'
      DOUBLE PRECISION Fx(*), Fy(*), Fz(*)
      DOUBLE PRECISION Fx2(*),Fy2(*), Fz2(*), En2, Vir2
      DOUBLE PRECISION Delt, dt, v2, Enkin

      INTEGER i, j, n

      v2 = 0.D0
      dt = Delt/n
      DO i = 1, NPART

         VX(i) = VX(i) + 0.5*Delt * (Fx(i)-Fx2(i))
         VY(i) = VY(i) + 0.5*Delt * (Fy(i)-Fy2(i))
         VZ(i) = VZ(i) + 0.5*Delt * (Fz(i)-Fz2(i))

         DO j = 1, n
            VX(i) = VX(i) + 0.5*dt * Fx2(i)
            VY(i) = VY(i) + 0.5*dt * Fy2(i)
            VZ(i) = VZ(i) + 0.5*dt * Fz2(i)

            X(i) = X(i) + dt * 2 * VX(i)
            Y(i) = Y(i) + dt * 2 * VY(i)
            Z(i) = Z(i) + dt * 2 * VZ(i)

            CALL FORCE(Fx2, Fy2, Fz2, En2, Vir2, nlist2, list2, 2)          !short

            VX(i) = VX(i) + 0.5*dt * Fx2(i)
            VY(i) = VY(i) + 0.5*dt * Fy2(i)
            VZ(i) = VZ(i) + 0.5*dt * Fz2(i)
         END DO

         CALL FORCE(Fx, Fy, Fz, En, Vir, nlist, list, 1)             !all

         VX(i) = VX(i) + 0.5*Delt * (Fx(i)-Fx2(i))
         VY(i) = VY(i) + 0.5*Delt * (Fy(i)-Fy2(i))
         VZ(i) = VZ(i) + 0.5*Delt * (Fz(i)-Fz2(i))

         v2 = v2 + (VX(i)**2 + VY(i)**2 + VZ(i)**2)
      END DO
      Enkin = v2/2
      RETURN
      END
