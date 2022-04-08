SUBROUTINE quantileMcpNETpath1 (gam,c,taux, lam2, maj, nobs, nvars, x, y, ju, &
& pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
        ! - - - arg types - - -
      DOUBLE PRECISION, PARAMETER :: big = 9.9E30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
      INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: ju (nvars)
      INTEGER :: m (pmax)
      INTEGER :: nbeta (nlam)
      INTEGER :: Nsort
      INTEGER :: NiSort
      INTEGER :: NvSort
      DOUBLE PRECISION :: sll
      DOUBLE PRECISION :: gam
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: c
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: SXij (nvars)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars)
    ! - - - local declarations - - -
      DOUBLE PRECISION :: d
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION :: u
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: al
      DOUBLE PRECISION :: alf
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: dl (nobs)
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
      INTEGER :: i
      INTEGER :: k
      INTEGER :: j
      INTEGER :: l
      INTEGER :: vrg
      INTEGER :: ctr
      INTEGER :: ierr
      INTEGER :: ni
      INTEGER :: me
      INTEGER, DIMENSION (:), ALLOCATABLE :: mm
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (b(0:nvars), STAT=jerr)
      ALLOCATE (oldbeta(0:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r(1:nobs), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
! - - - some initial setup - - -
      delta=c/max(taux,1-taux)
      r = y
      b = 0.0D0
      oldbeta = 0.0D0
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = Min (mnlam, nlam)
      maj = 2.0 * maj / delta
      IF (flmin < 1.0D0) THEN
         flmin = Max (mfl, flmin)
         alf = flmin ** (1.0D0/(nlam-1.0D0))
      END IF


! --------- lambda loop ----------------------------
      DO l = 1, nlam
         IF (flmin >= 1.0D0) THEN
            al = ulam (l)
         ELSE
            IF (l > 2) THEN
               al = al * alf
            ELSE IF (l == 1) THEN
               al = big
            ELSE IF (l == 2) THEN
               al = 0.0D0
               DO i = 1, nobs
                  IF (r(i) < (-c)) THEN
                     dl (i) = taux - 1.0D0
                  ELSE IF (r(i) < 0.0D0) THEN
                     dl (i) =  (1.0D0-taux) * r(i) / c
                  ELSE IF (r(i) < c) THEN
                     dl (i) = taux * r(i) / c
                  ELSE
                     dl (i) = taux
                  END IF
               END DO
               DO j = 1, nvars
                  IF (ju(j) /= 0) THEN
                     IF (pf(j) > 0.0D0) THEN
                        u = dot_product (dl, x(:, j))
                        al = Max (al, Abs(u)/pf(j))
                     END IF
                  END IF
               END DO
               al = al * alf / nobs
            END IF
         END IF
         ctr = 0
        ! --------- outer loop ----------------------------
         Nsort = 0
         DO
            Nsort = Nsort + 1
            IF (Nsort > 20) EXIT
            !!write(*,*),"boucle rouge-----------------------"
            oldbeta (0) = b (0)
            IF (ni > 0) oldbeta (m(1:ni)) = b (m(1:ni))
        ! --middle loop-------------------------------------
            Nvsort = 0
            DO
               Nvsort = Nvsort + 1
               IF (Nvsort > 20) EXIT
               !!write(*,*),"boucle vert_____________"
               npass = npass + 1
               dif = 0.0D0
               DO k = 1, nvars
                  IF (ju(k) /= 0) THEN
                     oldb = b (k)
                     u = 0.0D0
                    DO i = 1, nobs
                  	IF   (r(i) < (-c)) THEN
                     		dl (i) = taux-1.0D0
                  	ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0-taux) * r(i) / c
                  	ELSE IF (r(i) < c) THEN
                    		 dl (i) = taux * r(i) / c
                  	ELSE
                     		dl (i) = taux
                  	END IF
                        u = u + dl (i) * x (i, k)
                     END DO
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                    u = maj (k) * b (k) + u / nobs
                    v = al * pf (k)
                    IF (Abs(u) <= (gam * al * pf (k))) THEN
                         v = Abs (u) - v
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2 - 1.0/gam)
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE
                             b (k) = u/(maj(k) + pf2(k) * lam2)
                     END IF
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/delta)
                        r = r +  x (:, k) * d
                        IF (mm(k) == 0) THEN
                           ni = ni + 1
                           IF (ni > pmax) EXIT
                           mm (k) = ni
                           m (ni) = k !indicate which one is non-zero
                        END IF
                     END IF
                  END IF
               END DO

               IF (ni > pmax) EXIT

               d = 0.0D0
               DO i = 1, nobs
                 IF   (r(i) < (-c)) THEN
                     		dl (i) = taux-1.0D0
                 ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0-taux) * r(i) / c
                 ELSE IF (r(i) < c) THEN
                    		 dl (i) = taux * r(i) / c
                 ELSE
                     		dl (i) = taux
                 END IF
                        d = d + dl (i)
               END DO

               d = 0.5D0 * delta * d / nobs
               IF (d /= 0.0D0) THEN
                  b (0) = b (0) +  d
                  r = r - d
                  dif = Max (dif, 2.0*d**2/delta)
               END IF
               IF (dif < eps) EXIT
        ! --inner loop----------------------
               NiSort = 0
               DO
                  NiSort = NiSort + 1
                  IF (NiSort > 20) EXIT
                  npass = npass + 1
                  dif = 0.0D0
                  DO j = 1, ni
                     k = m (j)
                     oldb = b (k)
                     u = 0.0D0
                     DO i = 1, nobs
                 	IF (r(i) < -c) THEN
                     		dl (i) = taux-1.0D0
                 	ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0 - taux) * r(i) / c
                 	ELSE IF (r(i) < c) THEN
                    		dl (i) = taux * r(i) / c
                 	ELSE
                     		dl (i) = taux
                 	END IF
                        	u = u + dl (i) * x (i, k)
                     END DO
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                    u = maj (k) * b (k) + u / nobs
                    v = al * pf (k)
                    IF (Abs(u) <= (gam * al * pf (k))) THEN
                         v = Abs (u) - v
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2 - 1.0/gam)
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE
                             b (k) = u/(maj(k) + pf2(k) * lam2)
                     END IF
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/delta)
                        r = r + x (:, k) * d
                     END IF
                  END DO


               		d = 0.0D0
              		 DO i = 1, nobs
                	 IF   (r(i) < (-c)) THEN
                     		dl (i) = taux-1.0D0
                     ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0-taux) * r(i) / c
                 	ELSE IF (r(i) < c) THEN
                    		 dl (i) = taux * r(i) / c
                	 ELSE
                     		dl (i) = taux
                 	END IF
                        d = d + dl (i)
               		END DO

                   d = 0.5D0 * delta * d / nobs
                   IF (d /= 0.0D0) THEN
                     b (0) = b (0) +  d
                     r = r - d
                     dif = Max (dif, 2.0*d**2/delta)
                   END IF

                  IF (dif < eps) EXIT
                !!write(*,*),dif
               END DO
            END DO
            IF (ni > pmax) EXIT
        !--- this is the final check ------------------------
            vrg = 1
            IF ((b(0)-oldbeta(0))**2 >= eps) vrg = 0
            DO j = 1, ni
               IF ((b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
               jerr = - l
               RETURN
            END IF
          !!write(*,*),"fin boucle rouge=========================="
         END DO
    ! final update variable save results------------
         IF (ni > pmax) THEN
            jerr = - 10000 - l
            EXIT
         END IF
         IF (ni > 0) beta (1:ni, l) = b (m(1:ni))
         nbeta (l) = ni
         b0 (l) = b (0)
         alam (l) = al
         nalam = l
         IF (l < mnl) CYCLE
         IF (flmin >= 1.0D0) CYCLE
         me = count (beta(1:ni, l) /= 0.0D0)
         IF (me > dfmax) EXIT
      END DO
      DEALLOCATE (b, oldbeta, r, mm)
      RETURN
END SUBROUTINE quantileMcpNETpath1

! --------------------------------------------------
SUBROUTINE quantileMcpNET1 (gam, c,taux, lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: isd
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: jd (*)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: gam
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: c
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: ymean
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xmean(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xnorm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju)
      IF (jd(1) > 0) ju (jd(2:(jd(1)+1))) = 0

      IF (maxval(ju) <= 0) THEN
         jerr = 7777
         RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF

      pf = Max (0.0D0, pf)
      pf2 = Max (0.0D0, pf2)
      maj = 0.0D0
      CALL standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
      CALL quantileMcpNETpath1 (gam, c,taux, lam2, maj, nobs, nvars, x, y, ju, &
     & pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, &
     & ibeta, nbeta, alam, npass, jerr)
      IF (jerr > 0) RETURN! check error after calling function
! - - - organize beta afterward - - -
      DO l = 1, nalam
         nk = nbeta (l)
         IF (isd == 1) THEN
            DO j = 1, nk
               beta (j, l) = beta (j, l) / xnorm (ibeta(j))
            END DO
         END IF
          b0 (l) = b0 (l) - dot_product (beta(1:nk, l), &
          & xmean(ibeta(1:nk)))
      END DO
      DEALLOCATE (ju, xmean, xnorm, maj)
      RETURN
END SUBROUTINE quantileMcpNET1
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE quantileScadNETpath1 (gam,c,taux, lam2, maj, nobs, nvars, x, y, ju, &
& pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
        ! - - - arg types - - -
      DOUBLE PRECISION, PARAMETER :: big = 9.9E30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
      INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: Nsort
      INTEGER :: NiSort
      INTEGER :: NvSort
      INTEGER :: ju (nvars)
      INTEGER :: m (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: sll
      DOUBLE PRECISION :: gam
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: c
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: SXij (nvars)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars)
    ! - - - local declarations - - -
      DOUBLE PRECISION :: d
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION :: u
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: al
      DOUBLE PRECISION :: alf
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: dl (nobs)
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
      INTEGER :: i
      INTEGER :: k
      INTEGER :: j
      INTEGER :: l
      INTEGER :: vrg
      INTEGER :: ctr
      INTEGER :: ierr
      INTEGER :: ni
      INTEGER :: me
      INTEGER, DIMENSION (:), ALLOCATABLE :: mm
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (b(0:nvars), STAT=jerr)
      ALLOCATE (oldbeta(0:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r(1:nobs), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
! - - - some initial setup - - -
      delta=c/max(taux,1-taux)
      r = y
      b = 0.0D0
      oldbeta = 0.0D0
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = Min (mnlam, nlam)
      maj = 2.0 * maj / delta
      IF (flmin < 1.0D0) THEN
         flmin = Max (mfl, flmin)
         alf = flmin ** (1.0D0/(nlam-1.0D0))
      END IF


! --------- lambda loop ----------------------------
      DO l = 1, nlam
         IF (flmin >= 1.0D0) THEN
            al = ulam (l)
         ELSE
            IF (l > 2) THEN
               al = al * alf
            ELSE IF (l == 1) THEN
               al = big
            ELSE IF (l == 2) THEN
               al = 0.0D0
               DO i = 1, nobs
                  IF (r(i) < (-c)) THEN
                     dl (i) = taux - 1.0D0
                  ELSE IF (r(i) < 0.0D0) THEN
                     dl (i) =  (1.0D0-taux) * r(i) / c
                  ELSE IF (r(i) < c) THEN
                     dl (i) = taux * r(i) / c
                  ELSE
                     dl (i) = taux
                  END IF
               END DO
               DO j = 1, nvars
                  IF (ju(j) /= 0) THEN
                     IF (pf(j) > 0.0D0) THEN
                        u = dot_product (dl, x(:, j))
                        al = Max (al, Abs(u)/pf(j))
                     END IF
                  END IF
               END DO
               al = al * alf / nobs
            END IF
         END IF
         ctr = 0
        ! --------- outer loop ----------------------------
         Nsort = 0
         DO
            Nsort = Nsort + 1
            IF (Nsort > 20) EXIT
            oldbeta (0) = b (0)
            IF (ni > 0) oldbeta (m(1:ni)) = b (m(1:ni))
        ! --middle loop-------------------------------------
            Nvsort = 0
            DO
               Nvsort = Nvsort + 1
               IF (Nvsort > 20) EXIT
               npass = npass + 1
               dif = 0.0D0
               DO k = 1, nvars
                  IF (ju(k) /= 0) THEN
                     oldb = b (k)
                     u = 0.0D0
                    DO i = 1, nobs
                  	IF   (r(i) < (-c)) THEN
                     		dl (i) = taux-1.0D0
                  	ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0-taux) * r(i) / c
                  	ELSE IF (r(i) < c) THEN
                    		 dl (i) = taux * r(i) / c
                  	ELSE
                     		dl (i) = taux
                  	END IF
                        u = u + dl (i) * x (i, k)
                     END DO
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                    u = maj (k) * b (k) + u / nobs
                    v = al * pf (k)
                    IF (Abs(u) <= (al * pf (k))) THEN
                         v = Abs (u) - v
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2)
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE IF  (Abs(u) <= (gam*al * pf (k))) THEN
                         v = Abs (u) - gam*v/(gam-1)
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2-1/(gam-1))
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE
                         b (k) = u/(maj(k) + pf2(k) * lam2)
                     END IF
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/delta)
                        r = r +  x (:, k) * d
                        IF (mm(k) == 0) THEN
                           ni = ni + 1
                           IF (ni > pmax) EXIT
                           mm (k) = ni
                           m (ni) = k !indicate which one is non-zero
                        END IF
                     END IF
                  END IF
               END DO

               IF (ni > pmax) EXIT

               d = 0.0D0
               DO i = 1, nobs
                 IF   (r(i) < (-c)) THEN
                     		dl (i) = taux-1.0D0
                 ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0-taux) * r(i) / c
                 ELSE IF (r(i) < c) THEN
                    		 dl (i) = taux * r(i) / c
                 ELSE
                     		dl (i) = taux
                 END IF
                        d = d + dl (i)
               END DO

               d = 0.5D0 * delta * d / nobs
               IF (d /= 0.0D0) THEN
                  b (0) = b (0) +  d
                  r = r - d
                  dif = Max (dif, 2.0*d**2/delta)
               END IF

               IF (dif < eps) EXIT
        ! --inner loop----------------------
               NiSort = 0
               DO
                  NiSort = NiSort + 1
                  IF (NiSort > 20) EXIT
                  npass = npass + 1
                  dif = 0.0D0
                  DO j = 1, ni
                     k = m (j)
                     oldb = b (k)
                     u = 0.0D0
                     DO i = 1, nobs
                 	IF (r(i) < -c) THEN
                     		dl (i) = taux-1.0D0
                 	ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0 - taux) * r(i) / c
                 	ELSE IF (r(i) < c) THEN
                    		dl (i) = taux * r(i) / c
                 	ELSE
                     		dl (i) = taux
                 	END IF
                        	u = u + dl (i) * x (i, k)
                     END DO
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                    u = maj (k) * b (k) + u / nobs
                    v = al * pf (k)
                    IF (Abs(u) <= (al * pf (k))) THEN
                         v = Abs (u) - v
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2)
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE IF  (Abs(u) <= (gam*al * pf (k))) THEN
                         v = Abs (u) - gam*v/(gam-1)
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2-1/(gam-1))
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE
                         b (k) = u/(maj(k) + pf2(k) * lam2)
                     END IF
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/delta)
                        r = r + x (:, k) * d
                     END IF
                  END DO


               		d = 0.0D0
              		 DO i = 1, nobs
                	 IF   (r(i) < (-c)) THEN
                     		dl (i) = taux-1.0D0
                     ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0-taux) * r(i) / c
                 	ELSE IF (r(i) < c) THEN
                    		 dl (i) = taux * r(i) / c
                	 ELSE
                     		dl (i) = taux
                 	END IF
                        d = d + dl (i)
               		END DO

                   d = 0.5D0 * delta * d / nobs
                   IF (d /= 0.0D0) THEN
                     b (0) = b (0) +  d
                     r = r - d
                     dif = Max (dif, 2.0*d**2/delta)
                   END IF


                  IF (dif < eps) EXIT
               END DO
            END DO
            IF (ni > pmax) EXIT
        !--- this is the final check ------------------------
            vrg = 1
            IF ((b(0)-oldbeta(0))**2 >= eps) vrg = 0
            DO j = 1, ni
               IF ((b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
               jerr = - l
               RETURN
            END IF
         END DO
    ! final update variable save results------------
         IF (ni > pmax) THEN
            jerr = - 10000 - l
            EXIT
         END IF
         IF (ni > 0) beta (1:ni, l) = b (m(1:ni))
         nbeta (l) = ni
         b0 (l) = b (0)
         alam (l) = al
         nalam = l
         IF (l < mnl) CYCLE
         IF (flmin >= 1.0D0) CYCLE
         me = count (beta(1:ni, l) /= 0.0D0)
         IF (me > dfmax) EXIT
      END DO
      DEALLOCATE (b, oldbeta, r, mm)
      RETURN
END SUBROUTINE quantileScadNETpath1

! --------------------------------------------------
SUBROUTINE quantileScadNET1 (gam, c,taux, lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: isd
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: jd (*)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: gam
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: c
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: ymean
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xmean(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xnorm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju)
      IF (jd(1) > 0) ju (jd(2:(jd(1)+1))) = 0

      IF (maxval(ju) <= 0) THEN
         jerr = 7777
         RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF

      pf = Max (0.0D0, pf)
      pf2 = Max (0.0D0, pf2)
      maj = 0.0D0
      CALL standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
      CALL quantileScadNETpath1 (gam, c,taux, lam2, maj, nobs, nvars, x, y, ju, &
     & pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, &
     & ibeta, nbeta, alam, npass, jerr)
      IF (jerr > 0) RETURN! check error after calling function
! - - - organize beta afterward - - -
      DO l = 1, nalam
         nk = nbeta (l)
         IF (isd == 1) THEN
            DO j = 1, nk
               beta (j, l) = beta (j, l) / xnorm (ibeta(j))
            END DO
         END IF
          b0 (l) = b0 (l) - dot_product (beta(1:nk, l), &
          & xmean(ibeta(1:nk)))
      END DO
      DEALLOCATE (ju, xmean, xnorm, maj)
      RETURN
END SUBROUTINE quantileScadNET1
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE quantileScadNETpath2 (gam,kk,taux, lam2, maj, nobs, nvars, x, y, ju, &
& pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
        ! - - - arg types - - -
      DOUBLE PRECISION, PARAMETER :: big = 9.9E30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
      INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: ju (nvars)
      INTEGER :: m (pmax)
      INTEGER :: nbeta (nlam)
      INTEGER :: Nsort
      INTEGER :: NiSort
      INTEGER :: NvSort
      DOUBLE PRECISION :: gam
      DOUBLE PRECISION :: sll
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: kk
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: SXij (nvars)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars)
    ! - - - local declarations - - -
      DOUBLE PRECISION :: d
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: u
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: w
      DOUBLE PRECISION :: al
      DOUBLE PRECISION :: alf
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: dl (nobs)
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
      INTEGER :: i
      INTEGER :: k
      INTEGER :: j
      INTEGER :: l
      INTEGER :: vrg
      INTEGER :: ctr
      INTEGER :: ierr
      INTEGER :: ni
      INTEGER :: me
      INTEGER, DIMENSION (:), ALLOCATABLE :: mm
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (b(0:nvars), STAT=jerr)
      ALLOCATE (oldbeta(0:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r(1:nobs), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
! - - - some initial setup - - -
      r = y
      b = 0.0D0
      oldbeta = 0.0D0
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = Min (mnlam, nlam)
      maj = 2.0 * maj / kk
      IF (flmin < 1.0D0) THEN
         flmin = Max (mfl, flmin)
         alf = flmin ** (1.0D0/(nlam-1.0D0))
      END IF

! --------- lambda loop ----------------------------
      DO l = 1, nlam
         IF (flmin >= 1.0D0) THEN
            al = ulam (l)
         ELSE
            IF (l > 2) THEN
               al = al * alf
            ELSE IF (l == 1) THEN
               al = big
            ELSE IF (l == 2) THEN
               al = 0.0D0
               DO i = 1, nobs
                  IF (r(i) < (-taux*kk)) THEN
                     dl (i) = -taux
                  ELSE IF (r(i) < (1-taux)*kk) THEN
                     dl (i) =  r(i) / kk
                  ELSE
                     dl (i) = 1-taux
                  END IF
               END DO
               DO j = 1, nvars
                  IF (ju(j) /= 0) THEN
                     IF (pf(j) > 0.0D0) THEN
                        u = dot_product (dl, x(:, j))
                        al = Max (al, Abs(u)/pf(j))
                     END IF
                  END IF
               END DO
               al = al * alf / nobs
            END IF
         END IF
         ctr = 0
        ! --------- outer loop ----------------------------
         Nsort = 0
         DO
            Nsort = Nsort + 1
            IF (Nsort > 20) EXIT
            oldbeta (0) = b (0)
            IF (ni > 0) oldbeta (m(1:ni)) = b (m(1:ni))
        ! --middle loop-------------------------------------
            Nvsort = 0
            DO
               Nvsort = Nvsort + 1
               IF (Nvsort > 20) EXIT
               npass = npass + 1
               dif = 0.0D0
               DO k = 1, nvars
                  IF (ju(k) /= 0) THEN
                     oldb = b (k)
                     u = 0.0D0
                    DO i = 1, nobs
                 	IF (r(i) < (-taux*kk)) THEN
                     	dl (i) = -taux
                  	ELSE IF (r(i) < (1-taux)*kk) THEN
                     	dl (i) =  r(i) / kk
                  	ELSE
                     	dl (i) = 1-taux
                  	END IF
                        u = u + dl (i) * x (i, k)
                     END DO
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                    u = maj (k) * b (k) + u / nobs
                    v = al * pf (k)
                    IF (Abs(u) <= (al * pf (k))) THEN
                         v = Abs (u) - v
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2)
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE IF  (Abs(u) <= (gam*al * pf (k))) THEN
                         v = Abs (u) - gam*v/(gam-1)
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2-1/(gam-1))
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE
                         b (k) = u/(maj(k) + pf2(k) * lam2)
                     END IF
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/kk)
                        r = r +  x (:, k) * d
                        IF (mm(k) == 0) THEN
                           ni = ni + 1
                           IF (ni > pmax) EXIT
                           mm (k) = ni
                           m (ni) = k !indicate which one is non-zero
                        END IF
                     END IF
                  END IF
               END DO
               IF (ni > pmax) EXIT

                d = 0.0D0
				DO i = 1, nobs
                 	IF (r(i) < (-taux*kk)) THEN
                     	dl (i) = -taux
                  	ELSE IF (r(i) < (1-taux)*kk) THEN
                     	dl (i) =  r(i) / kk
                  	ELSE
                     	dl (i) = 1-taux
                  	END IF
                        d = d + dl (i)
                 END DO

                 d = 0.5D0 * kk * d / nobs
                 IF (d /= 0.0D0) THEN
                  b (0) = b (0) +  d
                  r = r - d
                  dif = Max (dif, 2.0*d**2/kk)
               END IF

               IF (dif < eps) EXIT
        ! --inner loop----------------------
               NiSort = 0
               DO
                  NiSort = NiSort + 1
                  IF (NiSort > 20) EXIT
                  npass = npass + 1
                  dif = 0.0D0
                  DO j = 1, ni
                     k = m (j)
                     oldb = b (k)
                     u = 0.0D0
                     DO i = 1, nobs
                  	 IF (r(i) < (-taux*kk)) THEN
                     		dl (i) = -taux
                  	 ELSE IF (r(i) < (1-taux)*kk) THEN
                     		dl (i) =  r(i) / kk
                  	 ELSE
                     		dl (i) = 1-taux
                  	 END IF
                        u = u + dl (i) * x (i, k)
                     END DO

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                    u = maj (k) * b (k) + u / nobs
                    v = al * pf (k)
                    IF (Abs(u) <= (al * pf (k))) THEN
                         v = Abs (u) - v
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2)
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE IF  (Abs(u) <= (gam*al * pf (k))) THEN
                         v = Abs (u) - gam*v/(gam-1)
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2-1/(gam-1))
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE
                         b (k) = u/(maj(k) + pf2(k) * lam2)
                     END IF
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/kk)
                        r = r + x (:, k) * d
                     END IF
                  END DO

               		d = 0.0D0
					DO i = 1, nobs
                 		IF (r(i) < (-taux*kk)) THEN
                     	dl (i) = -taux
                  		ELSE IF (r(i) < (1-taux)*kk) THEN
                     	dl (i) =  r(i) / kk
                  		ELSE
                     	dl (i) = 1-taux
                  		END IF
                        d = d + dl (i)
                	END DO

                   d = 0.5D0 * kk * d / nobs
                   IF (d /= 0.0D0) THEN
                     b (0) = b (0) +  d
                     r = r - d
                     dif = Max (dif, 2.0*d**2/kk)
                   END IF

                  IF (dif < eps) EXIT
               END DO
            END DO
            IF (ni > pmax) EXIT
        !--- this is the final check ------------------------
            vrg = 1
            IF ((b(0)-oldbeta(0))**2 >= eps) vrg = 0
            DO j = 1, ni
               IF ((b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
               jerr = - l
               RETURN
            END IF
         END DO
    ! final update variable save results------------
         IF (ni > pmax) THEN
            jerr = - 10000 - l
            EXIT
         END IF
         IF (ni > 0) beta (1:ni, l) = b (m(1:ni))
         nbeta (l) = ni
         b0 (l) = b (0)
         alam (l) = al
         nalam = l
         IF (l < mnl) CYCLE
         IF (flmin >= 1.0D0) CYCLE
         me = count (beta(1:ni, l) /= 0.0D0)
         IF (me > dfmax) EXIT
      END DO
      DEALLOCATE (b, oldbeta, r, mm)
      RETURN
END SUBROUTINE quantileScadNETpath2

! --------------------------------------------------
SUBROUTINE quantileScadNET2 (gam,kk,taux, lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: isd
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: jd (*)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: gam
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: kk
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: ymean
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xmean(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xnorm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju)
      IF (jd(1) > 0) ju (jd(2:(jd(1)+1))) = 0

      IF (maxval(ju) <= 0) THEN
         jerr = 7777
         RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF

      pf = Max (0.0D0, pf)
      pf2 = Max (0.0D0, pf2)
      maj = 0.0D0
      CALL standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
      CALL quantileScadNETpath2 (gam, kk,taux, lam2, maj, nobs, nvars, x, y, ju, &
     & pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, &
     & ibeta, nbeta, alam, npass, jerr)
      IF (jerr > 0) RETURN! check error after calling function
! - - - organize beta afterward - - -
      DO l = 1, nalam
         nk = nbeta (l)
         IF (isd == 1) THEN
            DO j = 1, nk
               beta (j, l) = beta (j, l)/ xnorm (ibeta(j))
            END DO
         END IF
          b0 (l) = b0 (l) - dot_product (beta(1:nk, l), &
          & xmean(ibeta(1:nk)))
      END DO
      DEALLOCATE (ju, xmean, xnorm, maj)
      RETURN
END SUBROUTINE quantileScadNET2
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------
SUBROUTINE quantileMcpNETpath2 (gam,kk,taux, lam2, maj, nobs, nvars, x, y, ju, &
& pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
        ! - - - arg types - - -
      DOUBLE PRECISION, PARAMETER :: big = 9.9E30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
      INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: ju (nvars)
      INTEGER :: m (pmax)
      INTEGER :: nbeta (nlam)
      INTEGER :: Nsort
      INTEGER :: NiSort
      INTEGER :: NvSort
      DOUBLE PRECISION :: gam
      DOUBLE PRECISION :: sll
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: kk
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: SXij (nvars)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars)
    ! - - - local declarations - - -
      DOUBLE PRECISION :: d
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: u
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: w
      DOUBLE PRECISION :: al
      DOUBLE PRECISION :: alf
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: dl (nobs)
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
      INTEGER :: i
      INTEGER :: k
      INTEGER :: j
      INTEGER :: l
      INTEGER :: vrg
      INTEGER :: ctr
      INTEGER :: ierr
      INTEGER :: ni
      INTEGER :: me
      INTEGER, DIMENSION (:), ALLOCATABLE :: mm
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (b(0:nvars), STAT=jerr)
      ALLOCATE (oldbeta(0:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r(1:nobs), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
! - - - some initial setup - - -
      r = y
      b = 0.0D0
      oldbeta = 0.0D0
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = Min (mnlam, nlam)
      maj = 2.0 * maj / kk
      IF (flmin < 1.0D0) THEN
         flmin = Max (mfl, flmin)
         alf = flmin ** (1.0D0/(nlam-1.0D0))
      END IF

! --------- lambda loop ----------------------------
      DO l = 1, nlam
         IF (flmin >= 1.0D0) THEN
            al = ulam (l)
         ELSE
            IF (l > 2) THEN
               al = al * alf
            ELSE IF (l == 1) THEN
               al = big
            ELSE IF (l == 2) THEN
               al = 0.0D0
               DO i = 1, nobs
                  IF (r(i) < (-taux*kk)) THEN
                     dl (i) = -taux
                  ELSE IF (r(i) < (1-taux)*kk) THEN
                     dl (i) =  r(i) / kk
                  ELSE
                     dl (i) = 1-taux
                  END IF
               END DO
               DO j = 1, nvars
                  IF (ju(j) /= 0) THEN
                     IF (pf(j) > 0.0D0) THEN
                        u = dot_product (dl, x(:, j))
                        al = Max (al, Abs(u)/pf(j))
                     END IF
                  END IF
               END DO
               al = al * alf / nobs
            END IF
         END IF
         ctr = 0
        ! --------- outer loop ----------------------------
         Nsort = 0
         DO
            Nsort = Nsort + 1
            IF (Nsort > 20) EXIT
            oldbeta (0) = b (0)
            IF (ni > 0) oldbeta (m(1:ni)) = b (m(1:ni))
        ! --middle loop-------------------------------------
            Nvsort = 0
            DO
               Nvsort = Nvsort + 1
               IF (Nvsort > 20) EXIT
               npass = npass + 1
               dif = 0.0D0
               DO k = 1, nvars
                  IF (ju(k) /= 0) THEN
                     oldb = b (k)
                     u = 0.0D0
                    DO i = 1, nobs
                 	IF (r(i) < (-taux*kk)) THEN
                     	dl (i) = -taux
                  	ELSE IF (r(i) < (1-taux)*kk) THEN
                     	dl (i) =  r(i) / kk
                  	ELSE
                     	dl (i) = 1-taux
                  	END IF
                        u = u + dl (i) * x (i, k)
                     END DO
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                    u = maj (k) * b (k) + u / nobs
                    v = al * pf (k)
                    IF (Abs(u) <= (gam * al * pf (k))) THEN
                         v = Abs (u) - v
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2 - 1.0/gam)
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE
                             b (k) = u/(maj(k) + pf2(k) * lam2)
                     END IF
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/kk)
                        r = r +  x (:, k) * d
                        IF (mm(k) == 0) THEN
                           ni = ni + 1
                           IF (ni > pmax) EXIT
                           mm (k) = ni
                           m (ni) = k !indicate which one is non-zero
                        END IF
                     END IF
                  END IF
               END DO
               IF (ni > pmax) EXIT

                d = 0.0D0
				DO i = 1, nobs
                 	IF (r(i) < (-taux*kk)) THEN
                     	dl (i) = -taux
                  	ELSE IF (r(i) < (1-taux)*kk) THEN
                     	dl (i) =  r(i) / kk
                  	ELSE
                     	dl (i) = 1-taux
                  	END IF
                        d = d + dl (i)
                 END DO

                 d = 0.5D0 * kk * d / nobs
                 IF (d /= 0.0D0) THEN
                  b (0) = b (0) +  d
                  r = r - d
                  dif = Max (dif, 2.0*d**2/kk)
               END IF

               IF (dif < eps) EXIT
        ! --inner loop----------------------
               NiSort = 0
               DO
                  NiSort = NiSort + 1
                  IF (NiSort > 20) EXIT
                  npass = npass + 1
                  dif = 0.0D0
                  DO j = 1, ni
                     k = m (j)
                     oldb = b (k)
                     u = 0.0D0
                     DO i = 1, nobs
                  	 IF (r(i) < (-taux*kk)) THEN
                     		dl (i) = -taux
                  	 ELSE IF (r(i) < (1-taux)*kk) THEN
                     		dl (i) =  r(i) / kk
                  	 ELSE
                     		dl (i) = 1-taux
                  	 END IF
                        u = u + dl (i) * x (i, k)
                     END DO

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                    u = maj (k) * b (k) + u / nobs
                    v = al * pf (k)
                    IF (Abs(u) <= (gam * al * pf (k))) THEN
                         v = Abs (u) - v
                         IF (v > 0.0D0) THEN
                     	     b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2 - 1.0/gam)
                         ELSE
                             b (k) = 0.0D0
                         END IF
                     ELSE
                             b (k) = u/(maj(k) + pf2(k) * lam2)
                     END IF
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/kk)
                        r = r + x (:, k) * d
                     END IF
                  END DO

               		d = 0.0D0
					DO i = 1, nobs
                 		IF (r(i) < (-taux*kk)) THEN
                     	dl (i) = -taux
                  		ELSE IF (r(i) < (1-taux)*kk) THEN
                     	dl (i) =  r(i) / kk
                  		ELSE
                     	dl (i) = 1-taux
                  		END IF
                        d = d + dl (i)
                	END DO

                   d = 0.5D0 * kk * d / nobs
                   IF (d /= 0.0D0) THEN
                     b (0) = b (0) +  d
                     r = r - d
                     dif = Max (dif, 2.0*d**2/kk)
                   END IF

                  IF (dif < eps) EXIT
               END DO
            END DO
            IF (ni > pmax) EXIT
        !--- this is the final check ------------------------
            vrg = 1
            IF ((b(0)-oldbeta(0))**2 >= eps) vrg = 0
            DO j = 1, ni
               IF ((b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
               jerr = - l
               RETURN
            END IF
         END DO
    ! final update variable save results------------
         IF (ni > pmax) THEN
            jerr = - 10000 - l
            EXIT
         END IF
         IF (ni > 0) beta (1:ni, l) = b (m(1:ni))
         nbeta (l) = ni
         b0 (l) = b (0)
         alam (l) = al
         nalam = l
         IF (l < mnl) CYCLE
         IF (flmin >= 1.0D0) CYCLE
         me = count (beta(1:ni, l) /= 0.0D0)
         IF (me > dfmax) EXIT
      END DO
      DEALLOCATE (b, oldbeta, r, mm)
      RETURN
END SUBROUTINE quantileMcpNETpath2

! --------------------------------------------------
SUBROUTINE quantileMcpNET2 (gam,kk,taux, lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: isd
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: jd (*)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: gam
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: kk
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: ymean
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xmean(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xnorm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju)
      IF (jd(1) > 0) ju (jd(2:(jd(1)+1))) = 0

      IF (maxval(ju) <= 0) THEN
         jerr = 7777
         RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF

      pf = Max (0.0D0, pf)
      pf2 = Max (0.0D0, pf2)
      maj = 0.0D0
      CALL standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
      CALL quantileMcpNETpath2 (gam, kk,taux, lam2, maj, nobs, nvars, x, y, ju, &
     & pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, &
     & ibeta, nbeta, alam, npass, jerr)
      IF (jerr > 0) RETURN! check error after calling function
! - - - organize beta afterward - - -
      DO l = 1, nalam
         nk = nbeta (l)
         IF (isd == 1) THEN
            DO j = 1, nk
               beta (j, l) = beta (j, l)/ xnorm (ibeta(j))
            END DO
         END IF
          b0 (l) = b0 (l) - dot_product (beta(1:nk, l), &
          & xmean(ibeta(1:nk)))
      END DO
      DEALLOCATE (ju, xmean, xnorm, maj)
      RETURN
END SUBROUTINE quantileMcpNET2
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------
SUBROUTINE standard(nobs,nvars,x,ju,isd,xmean,xnorm,maj)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    INTEGER::  nobs
    INTEGER::nvars
    INTEGER::isd
    INTEGER::ju(nvars)
    DOUBLE PRECISION::  x(nobs,nvars)
    DOUBLE PRECISION::xmean(nvars)
    DOUBLE PRECISION::xnorm(nvars)
    DOUBLE PRECISION::maj(nvars)
    ! - - - local declarations - - -
    INTEGER:: j
! - - - begin - - -
    DO j=1,nvars
        IF(ju(j)==1) THEN
            xmean(j)=sum(x(:,j))/nobs     !mean
            x(:,j)=x(:,j)-xmean(j)
            maj(j)=dot_product(x(:,j),x(:,j))/nobs
              IF(isd==1) THEN
                xnorm(j)=sqrt(maj(j))    !standard deviation
                x(:,j)=x(:,j)/xnorm(j)
                maj(j)=1.0D0
            ENDIF
        ENDIF
    ENDDO
END SUBROUTINE standard

! --------------------------------------------------
SUBROUTINE chkvars (nobs, nvars, x, ju)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: ju (nvars)
      DOUBLE PRECISION :: x (nobs, nvars)
    ! - - - local declarations - - -
      INTEGER :: i
      INTEGER :: j
      DOUBLE PRECISION :: t
! - - - begin - - -
      DO j = 1, nvars
         ju (j) = 0
         t = x (1, j)
         DO i = 2, nobs
            IF (x(i, j) /= t) THEN
               ju (j) = 1
               EXIT
            END IF
         END DO
      END DO
END SUBROUTINE chkvars
! --------------------------------------------------
SUBROUTINE quantilelassoNETpath1 (c,taux, lam2, maj, nobs, nvars, x, y, ju, &
& pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
        ! - - - arg types - - -
      DOUBLE PRECISION, PARAMETER :: big = 9.9E30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
      INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: ju (nvars)
      INTEGER :: m (pmax)
      INTEGER :: nbeta (nlam)
      INTEGER :: Nsort
      INTEGER :: NiSort
      INTEGER :: NvSort
      DOUBLE PRECISION :: sll
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: c
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: SXij (nvars)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars)
    ! - - - local declarations - - -
      DOUBLE PRECISION :: d
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION :: u
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: al
      DOUBLE PRECISION :: alf
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: dl (nobs)
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
      INTEGER :: i
      INTEGER :: k
      INTEGER :: j
      INTEGER :: l
      INTEGER :: vrg
      INTEGER :: ctr
      INTEGER :: ierr
      INTEGER :: ni
      INTEGER :: me
      INTEGER, DIMENSION (:), ALLOCATABLE :: mm
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (b(0:nvars), STAT=jerr)
      ALLOCATE (oldbeta(0:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r(1:nobs), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
! - - - some initial setup - - -
      delta=c/max(taux,1-taux)
      r = y
      b = 0.0D0
      oldbeta = 0.0D0
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = Min (mnlam, nlam)
      maj = 2.0 * maj / delta
      IF (flmin < 1.0D0) THEN
         flmin = Max (mfl, flmin)
         alf = flmin ** (1.0D0/(nlam-1.0D0))
      END IF


! --------- lambda loop ----------------------------
      DO l = 1, nlam
         IF (flmin >= 1.0D0) THEN
            al = ulam (l)
         ELSE
            IF (l > 2) THEN
               al = al * alf
            ELSE IF (l == 1) THEN
               al = big
            ELSE IF (l == 2) THEN
               al = 0.0D0
               DO i = 1, nobs
                  IF (r(i) < (-c)) THEN
                     dl (i) = taux - 1.0D0
                  ELSE IF (r(i) < 0.0D0) THEN
                     dl (i) =  (1.0D0-taux) * r(i) / c
                  ELSE IF (r(i) < c) THEN
                     dl (i) = taux * r(i) / c
                  ELSE
                     dl (i) = taux
                  END IF
               END DO
               DO j = 1, nvars
                  IF (ju(j) /= 0) THEN
                     IF (pf(j) > 0.0D0) THEN
                        u = dot_product (dl, x(:, j))
                        al = Max (al, Abs(u)/pf(j))
                     END IF
                  END IF
               END DO
               al = al * alf / nobs
            END IF
         END IF
         ctr = 0
        ! --------- outer loop ----------------------------
         Nsort = 0
         DO
            Nsort = Nsort + 1
            IF (Nsort > 20) EXIT
            oldbeta (0) = b (0)
            IF (ni > 0) oldbeta (m(1:ni)) = b (m(1:ni))
        ! --middle loop-------------------------------------
            Nvsort = 0
            DO
               Nvsort = Nvsort + 1
               IF (Nvsort > 20) EXIT
               npass = npass + 1
               dif = 0.0D0
               DO k = 1, nvars
                  IF (ju(k) /= 0) THEN
                     oldb = b (k)
                     u = 0.0D0
                    DO i = 1, nobs
                  	IF   (r(i) < (-c)) THEN
                     		dl (i) = taux-1.0D0
                  	ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0-taux) * r(i) / c
                  	ELSE IF (r(i) < c) THEN
                    		 dl (i) = taux * r(i) / c
                  	ELSE
                     		dl (i) = taux
                  	END IF
                        u = u + dl (i) * x (i, k)
                     END DO
                     u = maj (k) * b (k)  + u / nobs
                     v = al * pf (k)
                     v = Abs (u) - v
                     IF (v > 0.0D0) THEN
                     	b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2)
                     ELSE
                        b (k) = 0.0D0
                     END IF
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/delta)
                        r = r +  x (:, k) * d
                        IF (mm(k) == 0) THEN
                           ni = ni + 1
                           IF (ni > pmax) EXIT
                           mm (k) = ni
                           m (ni) = k !indicate which one is non-zero
                        END IF
                     END IF
                  END IF
               END DO

               IF (ni > pmax) EXIT

               d = 0.0D0
               DO i = 1, nobs
                 IF   (r(i) < (-c)) THEN
                     		dl (i) = taux-1.0D0
                 ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0-taux) * r(i) / c
                 ELSE IF (r(i) < c) THEN
                    		 dl (i) = taux * r(i) / c
                 ELSE
                     		dl (i) = taux
                 END IF
                        d = d + dl (i)
               END DO

               d = 0.5D0 * delta * d / nobs
               IF (d /= 0.0D0) THEN
                  b (0) = b (0) +  d
                  r = r - d
                  dif = Max (dif, 2.0*d**2/delta)
               END IF

               IF (dif < eps) EXIT
        ! --inner loop----------------------
               NiSort = 0
               DO
                  NiSort = NiSort + 1
                  IF (NiSort > 20) EXIT
                  npass = npass + 1
                  dif = 0.0D0
                  DO j = 1, ni
                     k = m (j)
                     oldb = b (k)
                     u = 0.0D0
                     DO i = 1, nobs
                 	IF (r(i) < -c) THEN
                     		dl (i) = taux-1.0D0
                 	ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0 - taux) * r(i) / c
                 	ELSE IF (r(i) < c) THEN
                    		dl (i) = taux * r(i) / c
                 	ELSE
                     		dl (i) = taux
                 	END IF
                        	u = u + dl (i) * x (i, k)
                     END DO
                     u = maj (k) * b (k)  + u / nobs
                     v = al * pf (k)
                     v = Abs (u) - v
                     IF (v > 0.0D0) THEN
                        b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2)
                     ELSE
                        b (k) = 0.0D0
                     END IF
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/delta)
                        r = r + x (:, k) * d
                     END IF
                  END DO


               		d = 0.0D0
              		 DO i = 1, nobs
                	 IF   (r(i) < (-c)) THEN
                     		dl (i) = taux-1.0D0
                     ELSE IF (r(i) < 0.0D0) THEN
                     		dl (i) =  (1.0D0-taux) * r(i) / c
                 	ELSE IF (r(i) < c) THEN
                    		 dl (i) = taux * r(i) / c
                	 ELSE
                     		dl (i) = taux
                 	END IF
                        d = d + dl (i)
               		END DO

                   d = 0.5D0 * delta * d / nobs
                   IF (d /= 0.0D0) THEN
                     b (0) = b (0) +  d
                     r = r - d
                     dif = Max (dif, 2.0*d**2/delta)
                   END IF


                  IF (dif < eps) EXIT
               END DO
            END DO
            IF (ni > pmax) EXIT
        !--- this is the final check ------------------------
            vrg = 1
            IF ((b(0)-oldbeta(0))**2 >= eps) vrg = 0
            DO j = 1, ni
               IF ((b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
               jerr = - l
               RETURN
            END IF
         END DO
    ! final update variable save results------------
         IF (ni > pmax) THEN
            jerr = - 10000 - l
            EXIT
         END IF
         IF (ni > 0) beta (1:ni, l) = b (m(1:ni))
         nbeta (l) = ni
         b0 (l) = b (0)
         alam (l) = al
         nalam = l
         IF (l < mnl) CYCLE
         IF (flmin >= 1.0D0) CYCLE
         me = count (beta(1:ni, l) /= 0.0D0)
         IF (me > dfmax) EXIT
      END DO
      DEALLOCATE (b, oldbeta, r, mm)
      RETURN
END SUBROUTINE quantilelassoNETpath1

! --------------------------------------------------
SUBROUTINE quantilelassoNET1 (c,taux, lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: isd
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: jd (*)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: c
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: ymean
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xmean(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xnorm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju)
      IF (jd(1) > 0) ju (jd(2:(jd(1)+1))) = 0

      IF (maxval(ju) <= 0) THEN
         jerr = 7777
         RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF

      pf = Max (0.0D0, pf)
      pf2 = Max (0.0D0, pf2)
      maj = 0.0D0
      CALL standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
      CALL quantilelassoNETpath1 (c,taux, lam2, maj, nobs, nvars, x, y, ju, &
     & pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, &
     & ibeta, nbeta, alam, npass, jerr)
      IF (jerr > 0) RETURN! check error after calling function
! - - - organize beta afterward - - -
      DO l = 1, nalam
         nk = nbeta (l)
         IF (isd == 1) THEN
            DO j = 1, nk
               beta (j, l) = beta (j, l) / xnorm (ibeta(j))
            END DO
         END IF
          b0 (l) = b0 (l) - dot_product (beta(1:nk, l), &
          & xmean(ibeta(1:nk)))
      END DO
      DEALLOCATE (ju, xmean, xnorm, maj)
      RETURN
END SUBROUTINE quantilelassoNET1

! --------------------------------------------------
SUBROUTINE quantilelassoNETpath2 (kk,taux, lam2, maj, nobs, nvars, x, y, ju, &
& pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
        ! - - - arg types - - -
      DOUBLE PRECISION, PARAMETER :: big = 9.9E30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
      INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: ju (nvars)
      INTEGER :: m (pmax)
      INTEGER :: nbeta (nlam)
      INTEGER :: Nsort
      INTEGER :: NiSort
      INTEGER :: NvSort
      DOUBLE PRECISION :: sll
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: kk
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: SXij (nvars)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars)
    ! - - - local declarations - - -
      DOUBLE PRECISION :: d
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: u
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: al
      DOUBLE PRECISION :: alf
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: dl (nobs)
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
      INTEGER :: i
      INTEGER :: k
      INTEGER :: j
      INTEGER :: l
      INTEGER :: vrg
      INTEGER :: ctr
      INTEGER :: ierr
      INTEGER :: ni
      INTEGER :: me
      INTEGER, DIMENSION (:), ALLOCATABLE :: mm
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (b(0:nvars), STAT=jerr)
      ALLOCATE (oldbeta(0:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r(1:nobs), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
! - - - some initial setup - - -
      r = y
      b = 0.0D0
      oldbeta = 0.0D0
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = Min (mnlam, nlam)
      maj = 2.0 * maj / kk
      IF (flmin < 1.0D0) THEN
         flmin = Max (mfl, flmin)
         alf = flmin ** (1.0D0/(nlam-1.0D0))
      END IF

! --------- lambda loop ----------------------------
      DO l = 1, nlam
         IF (flmin >= 1.0D0) THEN
            al = ulam (l)
         ELSE
            IF (l > 2) THEN
               al = al * alf
            ELSE IF (l == 1) THEN
               al = big
            ELSE IF (l == 2) THEN
               al = 0.0D0
               DO i = 1, nobs
                  IF (r(i) < (-taux*kk)) THEN
                     dl (i) = -taux
                  ELSE IF (r(i) < (1-taux)*kk) THEN
                     dl (i) =  r(i) / kk
                  ELSE
                     dl (i) = 1-taux
                  END IF
               END DO
               DO j = 1, nvars
                  IF (ju(j) /= 0) THEN
                     IF (pf(j) > 0.0D0) THEN
                        u = dot_product (dl, x(:, j))
                        al = Max (al, Abs(u)/pf(j))
                     END IF
                  END IF
               END DO
               al = al * alf / nobs
            END IF
         END IF
         ctr = 0
        ! --------- outer loop ----------------------------
         Nsort = 0
         DO
            Nsort = Nsort + 1
            IF (Nsort > 20) EXIT
            oldbeta (0) = b (0)
            IF (ni > 0) oldbeta (m(1:ni)) = b (m(1:ni))
        ! --middle loop-------------------------------------
            Nvsort = 0
            DO
               Nvsort = Nvsort + 1
               IF (Nvsort > 20) EXIT
               npass = npass + 1
               dif = 0.0D0
               DO k = 1, nvars
                  IF (ju(k) /= 0) THEN
                     oldb = b (k)
                     u = 0.0D0
                    DO i = 1, nobs
                 	IF (r(i) < (-taux*kk)) THEN
                     	dl (i) = -taux
                  	ELSE IF (r(i) < (1-taux)*kk) THEN
                     	dl (i) =  r(i) / kk
                  	ELSE
                     	dl (i) = 1-taux
                  	END IF
                        u = u + dl (i) * x (i, k)
                     END DO
                     u = maj (k) * b (k) + u / nobs
                     v = al * pf (k)
                     v = Abs (u) - v
                     IF (v > 0.0D0) THEN
                     	b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2)
                     ELSE
                        b (k) = 0.0D0
                     END IF
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/kk)
                        r = r +  x (:, k) * d
                        IF (mm(k) == 0) THEN
                           ni = ni + 1
                           IF (ni > pmax) EXIT
                           mm (k) = ni
                           m (ni) = k !indicate which one is non-zero
                        END IF
                     END IF
                  END IF
               END DO
               IF (ni > pmax) EXIT

                d = 0.0D0
				DO i = 1, nobs
                 	IF (r(i) < (-taux*kk)) THEN
                     	dl (i) = -taux
                  	ELSE IF (r(i) < (1-taux)*kk) THEN
                     	dl (i) =  r(i) / kk
                  	ELSE
                     	dl (i) = 1-taux
                  	END IF
                        d = d + dl (i)
                 END DO

                 d = 0.5D0 * kk * d / nobs
                 IF (d /= 0.0D0) THEN
                  b (0) = b (0) +  d
                  r = r - d
                  dif = Max (dif, 2.0*d**2/kk)
               END IF

               IF (dif < eps) EXIT
        ! --inner loop----------------------
               NiSort = 0
               DO
                  NiSort = NiSort + 1
                  IF (NiSort > 20) EXIT
                  npass = npass + 1
                  dif = 0.0D0
                  DO j = 1, ni
                     k = m (j)
                     oldb = b (k)
                     u = 0.0D0
                     DO i = 1, nobs
                  	IF (r(i) < (-taux*kk)) THEN
                     		dl (i) = -taux
                  	ELSE IF (r(i) < (1-taux)*kk) THEN
                     		dl (i) =  r(i) / kk
                  	ELSE
                     		dl (i) = 1-taux
                  	END IF
                        u = u + dl (i) * x (i, k)
                     END DO
                     u = maj (k) * b (k) + u / nobs
                     v = al * pf (k)
                     v = Abs (u) - v
                     IF (v > 0.0D0) THEN
                        b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2)
                     ELSE
                        b (k) = 0.0D0
                     END IF
                     d = oldb - b (k)
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/kk)
                        r = r + x (:, k) * d
                     END IF
                  END DO

               		d = 0.0D0
					DO i = 1, nobs
                 		IF (r(i) < (-taux*kk)) THEN
                     	dl (i) = -taux
                  		ELSE IF (r(i) < (1-taux)*kk) THEN
                     	dl (i) =  r(i) / kk
                  		ELSE
                     	dl (i) = 1-taux
                  		END IF
                        d = d + dl (i)
                	END DO

                   d = 0.5D0 * kk * d / nobs
                   IF (d /= 0.0D0) THEN
                     b (0) = b (0) +  d
                     r = r - d
                     dif = Max (dif, 2.0*d**2/kk)
                   END IF

                  IF (dif < eps) EXIT
               END DO
            END DO
            IF (ni > pmax) EXIT
        !--- this is the final check ------------------------
            vrg = 1
            IF ((b(0)-oldbeta(0))**2 >= eps) vrg = 0
            DO j = 1, ni
               IF ((b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
               jerr = - l
               RETURN
            END IF
         END DO
    ! final update variable save results------------
         IF (ni > pmax) THEN
            jerr = - 10000 - l
            EXIT
         END IF
         IF (ni > 0) beta (1:ni, l) = b (m(1:ni))
         nbeta (l) = ni
         b0 (l) = b (0)
         alam (l) = al
         nalam = l
         IF (l < mnl) CYCLE
         IF (flmin >= 1.0D0) CYCLE
         me = count (beta(1:ni, l) /= 0.0D0)
         IF (me > dfmax) EXIT
      END DO
      DEALLOCATE (b, oldbeta, r, mm)
      RETURN
END SUBROUTINE quantilelassoNETpath2

! --------------------------------------------------
SUBROUTINE quantilelassoNET2 (kk,taux, lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
& nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: isd
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: jd (*)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: kk
      DOUBLE PRECISION :: taux
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: ymean
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
! - - - begin - - -
! - - - allocate variables - - -
      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xmean(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xnorm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju)
      IF (jd(1) > 0) ju (jd(2:(jd(1)+1))) = 0

      IF (maxval(ju) <= 0) THEN
         jerr = 7777
         RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF

      pf = Max (0.0D0, pf)
      pf2 = Max (0.0D0, pf2)
      maj = 0.0D0
      CALL standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
      CALL quantilelassoNETpath2 (kk,taux, lam2, maj, nobs, nvars, x, y, ju, &
     & pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, &
     & ibeta, nbeta, alam, npass, jerr)
      IF (jerr > 0) RETURN! check error after calling function
! - - - organize beta afterward - - -
      DO l = 1, nalam
         nk = nbeta (l)
         IF (isd == 1) THEN
            DO j = 1, nk
               beta (j, l) = beta (j, l)/ xnorm (ibeta(j))
            END DO
         END IF
          b0 (l) = b0 (l) - dot_product (beta(1:nk, l), &
          & xmean(ibeta(1:nk)))
      END DO
      DEALLOCATE (ju, xmean, xnorm, maj)
      RETURN
END SUBROUTINE quantilelassoNET2
! -----------------------------------------------------------------------------------------------------------------------------

