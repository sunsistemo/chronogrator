      SUBROUTINE VLIST(nlist, list, rverlet, swiver)

c     swiver: switch to say which Verlet list we are using (short or long)

      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'system.inc'
      INCLUDE 'potential.inc'
      INCLUDE 'verlet.inc'

      INTEGER nlist(npmax), list(npmax, npmax)
      INTEGER i, j, swiver
      DOUBLE PRECISION xr, yr, zr, r, rverlet

      DO i = 1, NPART
         nlist(i) = 0
         IF (swiver.EQ.1) THEN
            XV(i) = X(i)
            YV(i) = Y(i)
            ZV(i) = Z(i)
         ELSE
            XV2(i) = X(i)
            YV2(i) = Y(i)
            ZV2(i) = Z(i)
         END IF
      END DO
      DO i = 1, NPART - 1
         DO j = i + 1, NPART
            xr = X(i) - X(j)
            yr = Y(i) - Y(j)
            zr = Z(i) - Z(j)

            IF (xr.GT.hbox) THEN
               xr = xr - box
            ELSE IF (xr.LT.-hbox) THEN
               xr = xr + box
            END IF
            IF (yr.GT.hbox) THEN
               yr = yr - box
            ELSE IF (yr.LT.-hbox) THEN
               yr = yr + box
            END IF
            IF (zr.GT.hbox) THEN
               zr = zr - box
            ELSE IF (zr.LT.-hbox) THEN
               zr = zr + box
            END IF

            r = SQRT(xr**2 + yr**2 + zr**2)
            IF (abs(r).LT.rverlet) THEN
               nlist(i) = nlist(i) + 1
               nlist(j) = nlist(j) + 1
               list(i, nlist(i)) = j
               list(j, nlist(j)) = i
            END IF
         END DO
      END DO
      RETURN
      END
