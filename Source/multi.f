      SUBROUTINE MULTI(Fx, Fy, Fz, Fx2, Fy2, Fz2, nlist, list,
     &                 nlist2, list2, Enkin, Enpot, Vir, Delt, n)

      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'veloc.inc'
      INCLUDE 'verlet.inc'
      DOUBLE PRECISION Fx(*), Fy(*), Fz(*)
      DOUBLE PRECISION Fx2(*),Fy2(*), Fz2(*), En2, Vir2
      DOUBLE PRECISION Delt, dt, v2, Enkin, Enpot, Vir, enpoti, viri
      INTEGER nlist(npmax), list(npmax, npmax)
      INTEGER nlist2(npmax), list2(npmax, npmax)

      INTEGER i, j, n

      v2 = 0.D0
      dt = Delt/n
      Enpot = 0.D0
      Enkin = 0.D0
      Vir = 0.D0

      DO i = 1, NPART
         VX(i) = VX(i) + 0.5*Delt * (Fx(i)-Fx2(i))
         VY(i) = VY(i) + 0.5*Delt * (Fy(i)-Fy2(i))
         VZ(i) = VZ(i) + 0.5*Delt * (Fz(i)-Fz2(i))
      END DO

      DO j = 1, n
         DO i = 1, NPART
            VX(i) = VX(i) + 0.5*dt * Fx2(i)
            VY(i) = VY(i) + 0.5*dt * Fy2(i)
            VZ(i) = VZ(i) + 0.5*dt * Fz2(i)

            X(i) = X(i) + dt * 2 * VX(i)
            Y(i) = Y(i) + dt * 2 * VY(i)
            Z(i) = Z(i) + dt * 2 * VZ(i)
            En2 = 0.D0
            Vir2 = 0.D0
            CALL f_particle(i, Fx2, Fy2, Fz2, En2, Vir2, nlist2, list2,
     &                      rv2, rdv2, 2)

            VX(i) = VX(i) + 0.5*dt * Fx2(i)
            VY(i) = VY(i) + 0.5*dt * Fy2(i)
            VZ(i) = VZ(i) + 0.5*dt * Fz2(i)
         END DO
      END DO

      DO i = 1, NPART
         enpoti = 0.D0
         viri = 0.D0
         CALL f_particle(i, Fx, Fy, Fz, enpoti, viri, nlist, list,
     &                   rv, rdv, 1)

         VX(i) = VX(i) + 0.5*Delt * (Fx(i)-Fx2(i))
         VY(i) = VY(i) + 0.5*Delt * (Fy(i)-Fy2(i))
         VZ(i) = VZ(i) + 0.5*Delt * (Fz(i)-Fz2(i))

         Vir = Vir + viri
         Enpot = Enpot + enpoti
         v2 = v2 + (VX(i)**2 + VY(i)**2 + VZ(i)**2)
      END DO
      Enkin = v2/2
      RETURN
      END
