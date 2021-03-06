**==md.spg  processed by SPAG 4.52O  at 15:46 on 28 Mar 1996
      PROGRAM MD
c________________________________________________________________________
c
c   Understanding Molecular Simulations: From Algorithms to Applications
c
c
c  We make no warranties, express or implied, that the programs contained
c  in this work are free of error, or that they will meet your requirements
c  for any particular application. They should not be relied on for solving
c  problems whose incorrect solution could results in injury, damage, or
c  loss of property. The authors and publishers disclaim all liability for
c  direct or consequential damages resulting from your use of the programs
c
c__________________________________________________________________________
c
c
c   Case Study 4:  Static properties of the Lennard-Jones fluid
c
c__________________________________________________________________________
c
      IMPLICIT NONE
      INTEGER nstep, nstep10, step, n
      LOGICAL scale, multi_on, drift
      DOUBLE PRECISION en, ent, vir, virt, enk, time, enpot, delt, tmax,
     &                 enkt, tequil, temprsq
      DOUBLE PRECISION enpot2, en2, vir2, enk2
c --common blocks declaration:
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'veloc.inc'
      INCLUDE 'potential.inc'
      INCLUDE 'system.inc'
      INCLUDE 'samp.inc'
      INCLUDE 'verlet.inc'
      INTEGER nlist(npmax), list(npmax, npmax)
      INTEGER nlist2(npmax), list2(npmax, npmax)
      DOUBLE PRECISION fx(NPMax), fy(NPMax), fz(NPMax)
      DOUBLE PRECISION fx2(NPMax), fy2(NPMax), fz2(NPMax)
      INTEGER E0_step
      DOUBLE PRECISION E0, E_drift

      E0_step = -1  ! steps since we recorded E_0 after equilibration
      E_drift = 0.D0

      WRITE (6, *) '**************** MC_NPT ***************'
c     ---initialize system
      CALL INIT(delt, tmax, tequil, temprsq, scale, multi_on, n, drift)
c     ---initialize Verlet lists
      CALL VLIST(nlist, list, rv, 1)
      IF (multi_on) THEN
         CALL VLIST(nlist2, list2, rv2, 2)
      END IF
c     --- Initial force calculation to kick-off multi time step integrator
      IF (multi_on) THEN
         CALL FORCE(Fx, Fy, Fz, Enpot, Vir, nlist, list, .true., 1)
         CALL FORCE(Fx2, Fy2, Fz2, Enpot2, Vir2, nlist2, list2, .true.,
     &        2)
      END IF
c     ---total energy of the system
      CALL TOTERG(en, vir, enk)
      WRITE (6, 99001) en - enk, enk, en, vir
      step = 0
      time = 0
      IF (SAMP1) CALL SAMPLE(0, step, en, vir, enk, delt)
      IF (SAMP2) CALL SAMPLE2(0, delt)
      nstep = INT(tmax/delt)
      nstep10 = INT(nstep/10)
      IF (nstep.EQ.0) nstep10 = 0
      DO WHILE (time.LT.tmax)
c     ---propagate all particles with one time-step and store new positions in a
c        Verlet-list
         IF (multi_on) THEN
            CALL MULTI(fx, fy, fz, fx2, fy2, fz2, nlist, list,
     &           nlist2, list2, enk, enpot, vir, delt, n)
            time = time + delt * n
         ELSE
            CALL FORCE(fx, fy, fz, enpot, vir, nlist, list, .false., 1)
            CALL SOLVE(fx, fy, fz, enk, delt)
            time = time + delt
         END IF
         en = enpot + enk
         step = step + 1
         IF (drift .and. E0_step.eq.-1 .and. time.gt.tequil) THEN
c     --- Set initial energy after equilibration
            CALL TOTERG(en2, vir2, enk2)
            E0 = en2
            E0_step = 0
         ELSE IF (drift .and. E0_step.gt.-1) THEN
c     --- Calculate energy drift
            CALL TOTERG(en2, vir2, enk2)
            E_drift = E_drift + abs((en2 - E0) / E0)
            E0_step = E0_step + 1
         END IF
         IF (time.LT.tequil) THEN
            IF (scale) THEN
               IF (MOD(step,20).EQ.0) CALL VELOCS(temprsq)
               IF (MOD(step,20).EQ.0) CALL VELOCS(temprsq)
            END IF
c           ---if system equilibrated sample averages:
         ELSE IF (MOD(step,NSAMP).EQ.0) THEN
            IF (SAMP1) CALL SAMPLE(1, step, en, vir, enk, delt)
            IF (SAMP2) CALL SAMPLE2(1, delt)
         END IF
         IF (MOD(step,nstep10).EQ.0) THEN
            WRITE (6, *) '======>> Done ', SNGL(time), ' out of ',
     &                   SNGL(tmax), en
c           ---write intermediate configuration to file
            CALL STORE(8)
         END IF
      END DO
      CALL TOTERG(ent, virt, enkt)
      IF (SAMP1) CALL SAMPLE(2, step, en, vir, enk, delt)
      IF (SAMP2) CALL SAMPLE2(2, delt)
      WRITE (6, 99002) ent - enkt, ent, virt
      CALL STORE(21)
c     --- Report average energy drift from the initial energy
      IF (drift) THEN
         E_drift = E_drift / (E0_step - 1)
         write (6, *) "Average Energy Drift: ", E_drift
      END IF
      STOP

99001 FORMAT (' Total pot. energy in. conf.       : ', f12.5, /,
     &        ' Total kinetic energy in. conf.    : ', f12.5, /,
     &        ' Total energy in. conf.            : ', f12.5, /,
     &        ' Total virial initial configuration: ', f12.5)
99002 FORMAT (' Potential erg. end of simulation  : ', f12.5, /,
     &     ' Total energy end of simulation       : ', f12.5, /,
     &     ' Total virial end of simulation       : ', f12.5)


      END
