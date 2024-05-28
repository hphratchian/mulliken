INCLUDE "mulliken_mod.f03"
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
      use mulliken_mod
!
!     Variable Declarations
!
      implicit none
      character(len=512)::matrixFilename
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      integer::numCmdLineArgs,MOanalysis,nAlpha,nBeta,nBasis,nBasisUse
      type(MQC_Variable)::nEalpha,nEbeta,nEtot
      type(MQC_Variable)::SMatrixAO,CMatrixAlpha,CMatrixBeta,  &
        PMatrixAlpha,PMatrixBeta,PMatrixTotal,basisFunctionTypeList
      type(MQC_Variable)::mullikenMatrixAO,mullikenGrossPopAO
      type(MQC_Variable)::tmpVec
!
!     Format Statements
!
 1000 Format(1x,'Enter Mulliken.')
 1010 Format(3x,'Matrix File: ',A)
 1020 Format(3x,'No MO specific analysis requested.',/)
 1022 Format(3x,'Analysis requested for MO number ',I5,'.',/)
 1100 Format(1x,'nAlpha=',I4,'  nBeta=',I4,'  nBasis=',I5,'  nBasisUse=',I5)
!
!
      write(IOut,1000)
!
!     Get the name of the matrix file and, optionally, get the number of the MO
!     for which an orbital Mulliken analysis is requested.
!
      numCmdLineArgs = command_argument_count()
      if(numCmdLineArgs.lt.1.or.numCmdLineArgs.gt.2)  &
        call mqc_error('Wrong number of command line arguments found. Use 1 or 2 arguments.')
      call get_command_argument(1,matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      MOanalysis = 0
      if(numCmdLineArgs.eq.2) call mqc_get_command_argument_integer(2,MOanalysis)
      if(MOanalysis.eq.0) then
        write(iOut,1020)
      else
        write(iOut,1022) MOanalysis
      endIf
!
!     Open the Gaussian matrix file and load the number of alpha electrons
!     (nAlpha), number of beta electrons (nBeta), and number of basis functions
!     (nBasis and nBasisUse).
!
      call GMatrixFile%load(matrixFilename)
      nAlpha    = GMatrixFile%getVal('nAlpha')
      nBeta     = GMatrixFile%getVal('nBeta')
      nBasis    = GMatrixFile%getVal('nBasis')
      nBasisUse = GMatrixFile%getVal('nBasisUse')
      basisFunctionTypeList = GMatrixFile%getBasisArray('basis type')
      write(IOut,1100) nAlpha,nBeta,nBasis,nBasisUse
!
!     Get the AO overlap and density matrices from the matrix file.
!
      call GMatrixFile%getArray('OVERLAP',mqcVarOut=SMatrixAO)
      call GMatrixFile%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=CMatrixAlpha)
      call GMatrixFile%getArray('ALPHA DENSITY MATRIX',mqcVarOut=PMatrixAlpha)
      if(GMatrixFile%isUnrestricted()) then
        call GMatrixFile%getArray('BETAMO COEFFICIENTS',mqcVarOut=CMatrixBeta)
        call GMatrixFile%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixBeta)
      else
        CMatrixBeta  = CMatrixAlpha
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
!     Carry out full Mulliken population analysis.
!
      call mulliken_analysis_ao(PMatrixTotal,SMatrixAO,mullikenMatrixAO,  &
        mullikenGrossPopAO)
      call mullikenMatrixAO%print(IOut,header='Mulliken Matrix in the AO Basis')
      call mullikenGrossPopAO%print(IOut,header='Gross AO Populations')
!
!     If requested, carry out a Mulliken population analysis on a specific
!     molecular orbital.
!

!hph
      if(MOanalysis.ne.0) then
        if(MOanalysis.gt.0) then
          tmpVec = CMatrixAlpha%column(MOanalysis)
        else
          tmpVec = CMatrixAlpha%column(-MOanalysis)
        endIf
        PMatrixTotal = outer_product(tmpVec,tmpVec)
        call PMatrixTotal%print(iOut,header='MO density matrix')
        call mulliken_analysis_ao(PMatrixTotal,SMatrixAO,mullikenMatrixAO,  &
          mullikenGrossPopAO,basisFunctionTypeList)
        call mullikenMatrixAO%print(IOut,header='Mulliken Matrix in the AO Basis for requested MO')
        call mullikenGrossPopAO%print(IOut,header='Gross AO Populations for requested MO')
        call mqc_print(sum(mullikenGrossPopAO),iOut,header='sum of gross pops',blankAtTop=.true.)
        mullikenGrossPopAO = mullikenGrossPopAO*mullikenGrossPopAO
        call mqc_print(sum(mullikenGrossPopAO),iOut,header='sum of gross pops SQUARED')
      endIf
!
!
!        procedure,public ::column        => MQC_Variable_MatrixGetColumn
!      function MQC_Variable_MatrixGetColumn(mqcVariable,iColumn)  &
!hph-


!
      end program mulliken
