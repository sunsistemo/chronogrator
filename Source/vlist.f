      SUBROUTINE VLIST(nlist, list)

      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'system.inc'
      INCLUDE 'potential.inc'
      INCLUDE 'verlet.inc'

      INTEGER nlist(npmax), list(npmax, npmax)
      INTEGER i, j
      DOUBLE PRECISION xr

      DO i = 1, NPART
         nlist(i) = 0
         XV(i) = X(i)
      END DO
      DO i = 1, NPART - 1
         DO j = i + 1, NPART
            xr = X(i) - X(j)
            IF (xr.GT.hbox) THEN
               xr = xr - box
            ELSE IF (xr.LT.-hbox) THEN
               xr = xr + box
            END IF
            IF (abs(xr).LT.rv) THEN
               nlist(i) = nlist(i) + 1
               nlist(j) = nlist(j) + 1
               list(i, nlist(i)) = j
               list(j, nlist(j)) = i
            END IF
         END DO
      END DO
      RETURN
      END
