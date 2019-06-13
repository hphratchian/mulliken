      program mulliken
!
!     This program carries out Mulliken population analysis from a Gaussian
!     unformatted matrix file.
!
!
!     H. P. Hratchian
!     Department of Chemistry & Chemical Biology
!     Center for Chemical Computation and Theory
!     University of California, Merced
!     hhratchian@ucmerced.edu
!
!
!
!     USE Connections
!
      use mqc_general
      use mqc_gaussian
      use mqc_algebra2
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
      character(len=512)::matrixFilename
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      integer::nAlpha,nBeta,nBasis,nBasisUse
      type(MQC_Variable)::nEalpha,nEbeta,nEtot
      type(MQC_Variable)::SMatrixAO,PMatrixAlpha,PMatrixBeta,PMatrixTotal
      type(MQC_Variable)::mullikenMatrixAO,mullikenGrossPopAO
!
!     Format Statements
!
 1000 Format(1x,'Enter Mulliken.')
 1010 Format(3x,'Matrix File: ',A,/)
 1100 Format(1x,'nAlpha=',I4,'  nBeta=',I4,'  nBasis=',I5,'  nBasisUse=',I5)
!
!
      write(IOut,1000)
!
!     Open the Gaussian matrix file and load the number of alpha electrons
!     (nAlpha), number of beta electrons (nBeta), and number of basis functions
!     (nBasis and nBasisUse).
!
      call get_command_argument(1,matrixFilename)
      call GMatrixFile%load(matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      nAlpha    = GMatrixFile%getVal('nAlpha')
      nBeta     = GMatrixFile%getVal('nBeta')
      nBasis    = GMatrixFile%getVal('nBasis')
      nBasisUse = GMatrixFile%getVal('nBasisUse')
      write(IOut,1100) nAlpha,nBeta,nBasis,nBasisUse
!
!     Get the AO overlap and density matrices from the matrix file.
!
      call GMatrixFile%getArray('OVERLAP',mqcVarOut=SMatrixAO)
      call GMatrixFile%getArray('ALPHA DENSITY MATRIX',mqcVarOut=PMatrixAlpha)
      if(GMatrixFile%isUnrestricted()) then
        call GMatrixFile%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixBeta)
      else
        PMatrixBeta  = PMatrixAlpha
      endIf
      PMatrixTotal = PMatrixAlpha+PMatrixBeta
      nEalpha = Contraction(PMatrixAlpha,SMatrixAO)
      nEbeta  = Contraction(PMatrixBeta,SMatrixAO)
      nEtot   = Contraction(PMatrixTotal,SMatrixAO)
      call nEalpha%print(IOut,'<P(Alpha)S>=')
      call nEbeta%print(IOut,'<P(Beta )S>=')
      call nEtot%print(IOut,'<P(Total)S>=')
!
!     Form the Mulliken matrix and compress it down to populationAO.
!
      mullikenMatrixAO = PMatrixTotal*SMatrixAO
      call mullikenMatrixAO%print(IOut,header='Mulliken Matrix in the AO Basis')
      
      mullikenGrossPopAO = MQC_Variable_Array_Sum(mullikenMatrixAO,1)
      call mullikenGrossPopAO%print(IOut,header='Gross AO Populations')

      end program mulliken
