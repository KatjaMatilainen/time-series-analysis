       PROGRAM DFT
C----------------------------------------------------------------------------
C  Author:   J. Lehar                                         Date: 13-JUL-87
C  Written for VAX FORTRAN 77   (May require some changes for other machines)
C  NOTE: This code was written for clarity, not for efficiency.
C        Some optimization may be desired for certain applications.
C----------------------------------------------------------------------------
C  Produces the dirty spectrum and the spectral window for a time series
C  contained in XFILE.  The time average and the data mean are automatically
C  removed (see paper).  The frequency sampling is determined by the user.
C              / freq. resolution = 1/T   ; T is the overall time interval
C   defaults: <  Max. frequency   = 1/2dt ; dt is the smallest time spacing
C              \ Points-per-beam  = 4
C  The dirty spectrum is written to DFILE, and the spectral window to WFILE.
C  NOTE: the filenames are limited to 10 characters.  This is so they can
C        be propagated to subsequent programs through the spec. file headers
C----------------------------------------------------------------------------
C  CALLS:     RDDATA       to read the time series data  \_ in READWRITE.FOR
C             WRSPEC       to write spectra              /
C             DFOUR        performs the discrete FT
C----------------------------------------------------------------------------
        PARAMETER (MAXN= 20000)           ! maximum # data samples
        PARAMETER (MAXM= 2000000)         ! maximum index for FT array

C  Declare arrays and some variables
        COMPLEX      D(0:MAXM),W(0:2*MAXM),DFOUR
        REAL         T(MAXN),X(MAXN),ONES(MAXN),F(0:2*MAXM)
        REAL         FREQ(MAXM),POWER(0:MAXM),WIN(MAXM),DFTAMP(0:MAXM)
        REAL         PMAX, FMAX
        CHARACTER    XFILE*20,DFILE*20,WFILE*20,HEADER*80
        DATA ONES/MAXN*1.0/           ! fill the ONES array with ones

C  Get execution parameters and filenames
         write(*,*) "INPUT FILES & PARAMETERS:"
       WRITE(*,100) ' Time-series data file: '
        READ (*,*) XFILE
       WRITE(*,100) ' Frequency resolution (0=default): '
        READ (*,*) FRES
       WRITE(*,100) ' Maximum frequency (0=default): '
        READ (*,*) FMAX
       WRITE(*,100) ' Points-per-beam (0=default): '
        READ (*,*) PPB
       WRITE (*,*) 
	 WRITE (*,*) ' Please, wait '
       WRITE (*,*) 

C  Read the times T(1:N) and the data X(1:N) from XFILE
        N= MAXN                               ! defines max. N for RDDATA
        CALL RDDATA(N,T,X,XFILE,HEADER)
        TMEAN= 0.
        XMEAN= 0.

C  If parameters are supplied as zero, set default values
        SMIN= 1.E20                          !\--(larger than expected SEPs)
        DO I=2,N                           ! \
           IF(T(I)-T(I-1).EQ.0) GOTO 898
           SEP= T(I) - T(I-1)                !  > find min. time separation
898        IF(SEP.LT.SMIN) SMIN= SEP         ! /
        ENDDO                                !/
        SMAX= T(N)-T(1)                      !> max. time separation
        IF(FRES.EQ.0.) FRES= 1./SMAX         ! frequency resolution
        IF(FMAX.EQ.0.) FMAX= 1./(2.*SMIN)    ! max. frequency
C        IF(PPB.EQ.0.) PPB= 4.                ! points-per-beam
        IF(PPB.EQ.0.) PPB= 1.                ! points-per-beam

C  set up the frequency array F(0:M)
        dF= FRES/PPB                   ! frequency increment
        M= INT(FMAX/dF)                ! maximum freq. element (for D)
C        type *, 'M=',M
        DO J=0,2*M
           F(J)= dF*J
        ENDDO

C  calculate the dirty spectrum D(0:M) and the spec. window W(0:2M)
C        TYPE *,'Computing the dirty spectrum...'
        DO J= 0,M
           D(J)= DFOUR(F(J),N,T,X)
        ENDDO
C        TYPE *,'Computing the spectral window...'
        DO J= 0,2*M
           W(J)= DFOUR(F(J),N,T,ONES)
        ENDDO


        OPEN(UNIT=1,FILE='psd.dat',STATUS='UNKNOWN')
        OPEN(UNIT=2,FILE='win.dat',STATUS='UNKNOWN')
        OPEN(UNIT=3,FILE='dft_amp.dat',STATUS='UNKNOWN')
C        J=0
        PMAXOLD=0.0
C        write the data to the file
        DO 22, J=0,M
C         DFTAMP(J)=SQRT(REAL(D(J))*REAL(D(J))+AIMAG(D(J))*AIMAG(D(J)))
           DFTAMP(J)=N*ABS(D(J))
           POWER(J)=2.*DFTAMP(J)**2/DFTAMP(0)
C           WIN(J)=SQRT(REAL(W(J))*REAL(W(J))+AIMAG(W(J))*AIMAG(W(J)))
           WIN(J)=N*ABS(W(J))
           IF (J.GT.0) THEN
              WRITE(1,'(3(1X,1P,E14.7))') F(J),POWER(J)
           ENDIF
           WRITE(2,'(3(1X,1P,E14.7))') F(J),WIN(J)
           WRITE(3,'(3(1X,1P,E14.7))') F(J),DFTAMP(J)
22      CONTINUE
        close(1)
        close(2)
        close(3)

C  exit the program
        CALL EXIT
100     Format(a,$)
        END




        COMPLEX FUNCTION DFOUR(FREQ,N,TIME,DATA)
C----------------------------------------------------------------------------
C  Author:   J. Lehar                                         Date: 13-JUL-87
C----------------------------------------------------------------------------
C  Returns the Fourier transform of the time series specified by TIME(1:N)
C  and DATA(1:N), evaluated at the frequency FREQ.
C  The form of the Fourier transform is taken from Bracewell (1965)
C
C                     1  .--. N           -i*2*PI*F*TIME(i)
C           DFT(F) = ---  >      DATA(i) e
C                     N  `--'i=1
C
C  The DFT is normalized to have the data mean at FREQ=0.
C----------------------------------------------------------------------------
        COMPLEX    SUM
        REAL       DATA(N),TIME(N)

C  initialize some variables
        SUM= (0.,0.)                    ! accumulation variable
        TWOPI= 8.*ATAN2(1.,1.)          ! calculate 2*pi

C  Evaluate FT at F...
        DO I=1,N
           PHASE= -TWOPI*FREQ*TIME(I)
           SUM= SUM + DATA(I)*CMPLX(COS(PHASE),SIN(PHASE))
        ENDDO

C  return with FT properly normalized
        DFOUR= SUM/N           ! ensures correct normalization
        RETURN
        END





        SUBROUTINE RDDATA(N,TIME,DATA,DFILE,HEADER)
C----------------------------------------------------------------------------
C  Author:  J. Lehar                                         Date:  13-JUL-87
C  Written for VAX FORTRAN 77   (May require some changes for other machines)
C  Source file:  READWRITE.FOR
C----------------------------------------------------------------------------
C  Reads in the times, data values for a time series from DFILE
C  In addition, a descriptive header is read in.  N (input) specifies the
C  maximum number of samples allowed and returns (output) the number of
C  samples found.  OPEN uses UNIT=1
C     ARGUMENTS:
C         N          (input) #samples allowed, (output) #samples found
C         TIME       (output) array for time samples
C         DATA       (output) array for data samples
C         DFILE      (input) data file name
C         HEADER     (output) descriptive file header
C----------------------------------------------------------------------------
        REAL TIME(N),DATA(N)
        CHARACTER DFILE*(*),HEADER*80

C  open the file & read in the header
        OPEN(UNIT=1,FILE=DFILE,STATUS='OLD')
C        READ(1,'(A80)') HEADER

C  read in the time series, each record has (Time,Data)
        I=1
C       step through the file until END-OF-FILE
c10         READ(1,'(2(1X,1P,E14.9))',END=11) TIME(I),DATA(I)
c10         READ(1,'(2(1X,E14.10))',END=11) TIME(I),DATA(I)
10          READ(1,*,END=11) TIME(I),DATA(I)
c           IF (TIME(I).EQ.TIME(I-1)) TYPE *, 'HERE', I, TIME(I)
c            IF (TIME(I).lt.9400.or.TIME(I).gt.9875) goto 10
C           TYPE *, TIME(I), DATA(I)
           I= I+1
           IF(I.GT.N)THEN
C              TYPE *,'RDDATA: data exceeds array size, truncating at N=',N
              GOTO 11
           ENDIF
           GOTO 10
11      CONTINUE
C        type *, 'here'
C  close the file & return the number of samples found
        CLOSE(UNIT=1)
        N= I-1
        RETURN
        END

