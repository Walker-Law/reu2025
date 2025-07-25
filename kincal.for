c     from https://userweb.jlab.org/~hibrahim/e01020/analysis/kincal/
c     Program: kincal.for
c     Purpose: Kinematics Calibration
c     Author : Hassan Ibrahim (2006)
c     Compile: g77 kincal.for -o kincal
c     Usage  : ./kincal

      implicit none

c     Declare Variables

      integer nx,ny,nzmax

      parameter (nx=6)
      parameter (ny=5)
      parameter (nzmax=20)       ! Set this >= number of runs

      double precision pi,mp

      parameter (pi=3.14159265359)
      parameter (mp=938.272)

      integer i,j,k,l,n,nz,lun,dof
      integer run(nzmax)

      double precision chisq,chisqpdof,chisqmin,chisqminpdof
      double precision e0(nzmax),pe0(nzmax),pp0(nzmax)
      double precision the0(nzmax),thp0(nzmax)
      double precision e,pe,pp,ep
      double precision sthe,cthe,sthp,cthp
      double precision sphe,cphe,sphp,cphp
      double precision sigx(nx),delx(nx),sigdelx(nx)
      double precision dely(ny,nzmax),delysum(ny),delyave(ny)
      double precision delynew(ny,nzmax)
      double precision delymin(ny,nzmax),delyminsum(ny),delyminave(ny)
      double precision sigy2(ny,nzmax)
      double precision dy_dx(nx,ny,nzmax)
      double precision mat(nx,nx),mat_inv(nx,nx)
      double precision vec(nx),vec_out(nx)

      character*80 dummy

c     Read initial offsets and nominal uncertainites

      open(unit=2,name='kincal.inp',status='old',form='formatted')

      read(2,*) dummy
      read(2,*) dummy
      read(2,*) dummy

      read(2,*) (sigx(i), i = 1, nx)
 
      read(2,*) dummy
      read(2,*) dummy
      read(2,*) dummy

      k = 1

 100  read(2,*,end=999) run(k),e0(k),pe0(k),pp0(k),the0(k),thp0(k),
     #        (dely(j,k), j = 1, ny)

      k = k + 1

      goto 100

 999  close(unit=2)

      nz = k - 1

c     Initialization

      chisq = 0.0

      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
               dy_dx(i,j,k) = 0.0
            enddo
         enddo
      enddo

      do j = 1, ny
         do k = 1, nz
            sigy2(j,k) = 0.0
            delynew(j,k) = 0.0
         enddo
      enddo

      do i = 1, nx
         vec(i) = 0.0
      enddo

      do i = 1, nx
         do l = 1, nx
            mat(i,l) = 0.0
         enddo
      enddo

      do j = 1, ny
         delysum(j) =0.0
      enddo

c     Start the main loop over kinematic settings

      do k = 1, nz

         e = e0(k)
         pe = pe0(k)
         pp = pp0(k)
         sthe = sin(the0(k) * pi / 180.0)
         cthe = cos(the0(k) * pi / 180.0)
         sthp = sin(thp0(k) * pi / 180.0)
         cthp = cos(thp0(k) * pi / 180.0)
         sphe = 0.0
         cphe = 1.0
         sphp = 0.0
         cphp = 1.0         
         ep = sqrt(pp**2 + mp**2)

         dy_dx(1,1,k) = -e                      ! dw_ddele  
         dy_dx(3,1,k) = -e*pe/mp*sthe*cphe      ! dw_dthe  
         dy_dx(5,1,k) = -e*pe/mp*cthe*sphe      ! dw_dphe
                 
         dy_dx(1,2,k) = -pe                     ! dem_ddele 
         dy_dx(2,2,k) = -pp**2/ep               ! dem_ddelp 
                 
         dy_dx(1,3,k) = -pe*sthe*cphe           ! dpmx_ddele
         dy_dx(2,3,k) = -pp*sthp*cphp           ! dpmx_ddelp
         dy_dx(3,3,k) = -pe*cthe*cphe           ! dpmx_dthe 
         dy_dx(4,3,k) = -pp*cthp*cphp           ! dpmx_dthp 
         dy_dx(5,3,k) = +pe*sthe*sphe           ! dpmx_dphe 
         dy_dx(6,3,k) = +pp*sthp*sphp           ! dpmx_dphp 
                 
         dy_dx(1,4,k) = -pe*sphe                ! dpmy_ddele
         dy_dx(2,4,k) = -pp*sphp                ! dpmy_ddelp
         dy_dx(5,4,k) = -pe*cphe                ! dpmy_dphe
         dy_dx(6,4,k) = -pp*cphp                ! dpmy_dphp
                 
         dy_dx(1,5,k) = -pe*cthe*cphe           ! dpmz_ddele
         dy_dx(2,5,k) = -pp*cthp*cphp           ! dpmz_ddelp
         dy_dx(3,5,k) = +pe*sthe*cphe           ! dpmz_dthe
         dy_dx(4,5,k) = +pp*sthp*cphp           ! dpmz_dthp
         dy_dx(5,5,k) = +pe*cthe*sphe           ! dpmz_dphe
         dy_dx(6,5,k) = +pp*cthp*sphp           ! dpmz_dphp   
     
         do i = 1, nx
            do j = 1, ny
               sigy2(j,k) = sigy2(j,k) + dy_dx(i,j,k)**2 * sigx(i)**2
            enddo
         enddo

         do i = 1, nx
            do j = 1, ny

               vec(i) = vec(i) + dely(j,k) * dy_dx(i,j,k) / sigy2(j,k)

               do l = 1, nx

                  mat(i,l) = mat(i,l) + dy_dx(i,j,k) * dy_dx(l,j,k) / 
     #                 sigy2(j,k)

               enddo
            enddo
         enddo

      enddo

c     End of the main loop over kinematic settings

      dof = (ny - 1) * (nz - 1)

c     Initial Chi-Square

      do j = 1, ny
         do k = 1, nz

            chisq = chisq +(dely(j,k))**2/sigy2(j,k)

         enddo
      enddo

      chisqpdof = chisq / dof

      do j = 1, ny
         do k = 1, nz
            delysum(j) = delysum(j) + dely(j,k)
         enddo
         delyave(j) = delysum(j) / nz
      enddo

c     Chi-Square Minimization

      call matinv(mat,mat_inv,nx)
      call mat_mult(nx,nx,1,mat_inv,vec,vec_out)

      do i = 1, nx
         delx(i) = vec_out(i)
         sigdelx(i) = sqrt(mat_inv(i,i))
      enddo

c     Minimized Chi-Square
      
      do j = 1, ny
         do k = 1, nz

            do i = 1, nx
               delynew(j,k) = delynew(j,k) + dy_dx(i,j,k) * delx(i)
            enddo

            delymin(j,k) = dely(j,k)-delynew(j,k)

            chisqmin = chisqmin +delymin(j,k)**2/sigy2(j,k)
            
         enddo
      enddo

      chisqminpdof = chisqmin / dof

      do j = 1, ny
         do k = 1, nz
            delyminsum(j) = delyminsum(j) + delymin(j,k)
         enddo
         delyminave(j) = delyminsum(j) / nz
      enddo

c     Write out the results

      do n = 1, 2

         if (n .eq. 1) then
            lun = 1
            open(unit=lun,name='kincal.out',status='unknown')
         else
            lun = 6
         endif
         write(lun,*)
         write(lun,*) "Kinematics Calibration"
         write(lun,*) "----------------------"
         write(lun,*)
         write(lun,*) "Nominal Spectrometer Uncertainties:"
         write(lun,*)
         write(lun,10) "DelE","DelP","TheE (rad)","TheP (rad)",
     #        "PhiE (rad)","PhiP (rad)"
         write(lun,*)
         write(lun,20) (sigx(i), i = 1, nx)
         write(lun,*)
         write(lun,*) "Number of Runs =",nz
         write(lun,*)
         write(lun,*) "Initial Kinematical Offsets:"
         write(lun,*)
         write(lun,30) "Run", "W (MeV)","Em (MeV)","Pmx (MeV)",
     #        "Pmy (MeV)", "Pmz (MeV)"
         write(lun,*)

         do k = 1, nz
            write(lun,40) run(k), (dely(j,k), j = 1, ny)
         enddo

         write(lun,*)
         write(lun,50) "AVG", (delyave(j),j=1,ny)
         write(lun,*)
         write(lun,*) "Initial Chi-Square per Degree of Freedom =",
     #        chisqpdof
         write(lun,*)
         write(lun,*) "------------------------------------------------"
         write(lun,*)
         write(lun,*) "Fit Results:"
         write(lun,*)
         write(lun,*) "Spectrometer Offsets:"
         write(lun,*)
         write(lun,10) "DelE","DelP","TheE (rad)","TheP (rad)",
     #        "PhiE (rad)","PhiP (rad)"
         write(lun,20) (delx(i), i = 1, nx)
         write(lun,*)


         write(lun,*) "Spectrometer Uncertainties:"
         write(lun,*)
         write(lun,20) (sigdelx(i), i = 1, nx)
         write(lun,*)
         write(lun,*) "------------------------------------------------"
         write(lun,*)
         write(lun,*) "Minimized Kinematical Offsets:"
         write(lun,*)
         write(lun,30) "Run", "W (MeV)","Em (MeV)","Pmx (MeV)",
     #        "Pmy (MeV)", "Pmz (MeV)"
         write(lun,*)

         do k = 1, nz
            write(lun,40) run(k), (delymin(j,k),j=1,ny)
         enddo

         write(lun,*)
         write(lun,50) "AVG", (delyminave(j),j=1,ny)
         write(lun,*)
         write(lun,*) "Minimized Chi-Square per Degree of Freedom =",
     #        chisqminpdof
         write(lun,*)

         close(unit=lun)

 10      format(6(a10,2x))
 20      format(6(f10.6,2x))
 30      format(a4,6(2x,a10))
 40      format(i4,5(2x,f10.6))
 50      format(a4,5(2x,f10.6))

      enddo

c     
      stop
      end
C     
C     --------------------------------------------------------------------
C     --------------------------------------------------------------------
C     SUBROUTINE MATINV
C     
C     Purpose:
C     Invert a symmetric matrix
C     
C     Usage:
C     CALL MATINV( ARRAY, ARRAY_INV, NORDER )
C     
C     Description of parameters:
C     ARRAY     - input matrix
C     ARRAY_INV - inverse matrix
C     NORDER    - degree of matrix
C     
C     Subroutines and function subprograms required:
C     none
C     
C     Based on routine of Bevington.
C     
C     --------------------------------------------------------------------
C     
      SUBROUTINE MATINV (ARRAY, ARRAY_INV, NORDER  )
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ARRAY(NORDER,NORDER),ARRAY_INV(NORDER,NORDER)
      INTEGER IK(20), JK(20)
C     
      DO I=1,NORDER
         DO J=1,NORDER
            ARRAY_INV(I,J) = ARRAY(I,J)
         ENDDO
      ENDDO
C     
      DO K=1,NORDER
C     
C     FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX
C     
         AMAX = 0.D0
 21      DO I=K,NORDER
            DO J=K, NORDER
               IF(  ABS(AMAX) -  ABS(ARRAY_INV(I,J) ) ) 24,24,30
 24            AMAX = ARRAY_INV( I,J )
               IK(K) = I
               JK(K) = J
 30         ENDDO
         ENDDO
C     
C     INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAY_INV(K,K)
C     
         IF( AMAX ) 41,140, 41
 41      I=IK(K)
         IF( I-K) 21,51,43
C     
 43      DO J=1,NORDER
            SAVE = ARRAY_INV( K,J)
            ARRAY_INV(K,J) = ARRAY_INV(I,J)
            ARRAY_INV(I,J) = -SAVE
         ENDDO
C     
 51      J = JK(K)
         IF(J-K) 21,61,53
C     
 53      DO I=1,NORDER
            SAVE = ARRAY_INV(I,K)
            ARRAY_INV(I,K) = ARRAY_INV(I,J)
            ARRAY_INV(I,J) = -SAVE
         ENDDO
C     
C     ACCUMULATE ELEMENTS OF INVERSE MATRIX
C     
 61      DO I=1,NORDER
            IF(I-K) 63,70,63
 63         ARRAY_INV(I,K) = -ARRAY_INV(I,K)/AMAX
 70      ENDDO
C     
 71      DO I=1,NORDER
            DO J=1,NORDER
               IF( I-K) 74,80,74
 74            IF( J-K) 75,80,75
 75            ARRAY_INV(I,J) = ARRAY_INV(I,J)
     #              + ARRAY_INV(I,K)*ARRAY_INV(K,J)
 80         ENDDO
         ENDDO
C     
 81      DO J=1,NORDER
            IF(J-K) 83,90,83
 83         ARRAY_INV(K,J) = ARRAY_INV(K,J)/AMAX
 90      ENDDO
C     
         ARRAY_INV(K,K) = 1./AMAX
      ENDDO
C     
C     RESTORE ORDERING OF MATRIX
C     
 101  DO L=1, NORDER
         K = NORDER-L+1
         J = IK(K)
         IF(J-K) 111,111,105
C     
 105     DO I=1,NORDER
            SAVE = ARRAY_INV(I,K)
            ARRAY_INV(I,K) = -ARRAY_INV(I,J)
            ARRAY_INV(I,J) = SAVE
         ENDDO
C     
 111     I = JK(K)
         IF(I-K) 130,130,113
C     
 113     DO J=1,NORDER
            SAVE = ARRAY_INV(K,J)
            ARRAY_INV(K,J) = -ARRAY_INV(I,J)
            ARRAY_INV(I,J) = SAVE
         ENDDO
 130  ENDDO
 140  RETURN
      END
C     
C     ----------------------------------------------------------------------
C     
C     SUBROUTINE MAT_MULT
C     AUTHOR:   M. Nozar
C     DATE:     19-JUL-1991
C     PURPOSE:
C     Computes product of any two matrices:
C     matprod(mxp) = matrix1(mxn) * matrix2(nxp)         
C     -----------------------------------------------------------------------
C     
      SUBROUTINE MAT_MULT(m,n,p,MATRIX1,MATRIX2,MATPROD)
      IMPLICIT NONE
C     
      INTEGER            m,n,p,I,J,K  
      DIMENSION          MATRIX1(m,n),MATRIX2(n,p),MATPROD(m,p)
      DOUBLE PRECISION   MATRIX1, MATRIX2, MATPROD
C     
C     ----------------------------------------------------------------------
C     Multiply two matrices
C     ----------------------------------------------------------------------
C     
      DO I = 1,m
         DO J = 1,p
            MATPROD(I,J) = 0.D0
            DO K = 1,n
               MATPROD(I,J) = MATPROD(I,J)+MATRIX1(I,K)*MATRIX2(K,J) 
            ENDDO         
         ENDDO
      ENDDO
C     
      RETURN
      END