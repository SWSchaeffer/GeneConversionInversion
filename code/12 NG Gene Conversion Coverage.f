      program nggcco
c
c-----NG Gene Conversion Coverage
c-----Estimates gene conversion tract coverage in 14 regions
c-----in coding sequence exons (CDS),introns, and intergenic regions
c-----file:///C:/Users/sws4/Documents/Computer Code/Fortran/NG Gene Conversion Permutation Test.f
c-----Written by Stephen W. Schaeffer
c-----Date: 8 May 2024
c
c-----Source: NG Gene Conversion Coverage.f
c
      integer ib,ic,id,ie,ig
      integer chr(19787792),rchr(19787792),igchr(19787792)
      integer nbeg(35000),nend(35000),nib,nie
      integer nr,rbeg(14),rend(14),rcov(14),rcov2(14)
      integer gcov(2,14),ggc,ggcs
c      integer gcovi(2,2,14)
      integer imin,imax
      real bp,acov(14),acov2(14)
      real bcov(2,14),gsit(2,14)
c      real gsiti(2,2,14),bcovi(2,2,14)
      integer ntr,ntb(11212),nte(11212),ntex(11212),ts(11212)
      character*1 t
      character*2 arr
      character*3 a
c      character*9 rname(11213)
      character*11 gnamt,gname(11212)
c
c-----nr  - number of syntenic regions
c-----arr - two letter code for the gene arrangement, user supplied
c
     	t=char(9)
      write(*,'("Input the inversion.")')
      read(*,*) arr
      nr=14
c
c-----Input data from the XX_MinMax.tsv file
c
c-----The file has coordinates of the 14 syntenic block regions
c-----rbeg(i) - beginning coordinate for the ith syntenic region
c-----rend(i) - end coordinate for the ith syntenic region
c-----rchr    - array where each element is an index with region number for each site
c
      open(unit=1,file=arr//'_MinMax.tsv')
      read(1,'(1x)')
      do i=1,nr
      read(1,*) rbeg(i),rend(i)
      end do
      close(unit=1)
      do i=1,nr
      do j=rbeg(i),rend(i)
      rchr(j)=i
      end do
      end do
c
c-----chr() - array of nucleotide positions on the third chromosome map
c-----Initialize chr()
c
      chr=0
c
c-----Open file with gene conversion tracts and input the data
c
c-----ngc      - number of gene conversion tracts
c-----nbeg(i)  - coordinate of the first nucleotide of the ith gene conversion tract
c-----nend(i)  - coordinate of the last  nucleotide in the ith gene conversion tract
c
      open(unit=1,file='GC_'//arr//'_Tracts.tsv')
      read(1,'(1x)')
c
c-----Read in the ngc coordinate pairs nbeg, nend until End Of File
c
      ngc=0      
50    ngc=ngc+1
      read(1,*,end=100) a,ib,ic,id,ie,ig,nbeg(ngc),nend(ngc)
c
c-----Increment elements in chr() using the gene conversion tract coodinates
c
c-----A value of 0 indicates no coverage, while a value > 0 indicates
c-----how many gene conversion tracts cover the nucleotide
c      
      do i=nbeg(ngc),nend(ngc)
      chr(i)=chr(i)+1
      end do
      goto 50      
100   close(unit=1)
      ngc=ngc-1
c
c-----Estimate gene conversion tract depth in genes
c
c-----ntr      - Number of exons across all genes
c-----gname(i) - Name of the ith transcript exon
c-----ntex(i)  - Number of the ith transcript exon
c-----ntb(i)   - First nucleotide of the ith transcript exon
c-----nte(i)   - Last nucleotide of the ith transcript exon
c-----ts(i)    - Status of the ith transcript exon [1=non-outlier; 2=outlier]
c      
      open(unit=1,file='GC_'//arr//'_Transcripts_Exons_List.tsv')
      ntr=0
150   ntr=ntr+1
      read(1,*,end=200) gname(ntr),ntex(ntr),ntb(ntr),nte(ntr),ts(ntr)
      goto 150
200   close(unit=1)
      ntr=ntr-1
c
c-----Estimate mean and variance of gene conversion coverage in each region
c-----rcov(i)  - sum of gene conversion tracts per nucleotide in the ith region
c-----rcov2(i) - sum of squares of gene conversion tracts per nucleotide in the ith region
c-----acov(i)  - mean gene conversion per nucleotide in the ith region
c-----acov2(i) - variance of gene conversion per nucleotide in the ith region
c
      rcov=0
      rcov2=0   
      do i=1,nr
      do j=rbeg(i),rend(i)
      rcov(i)=rcov(i)+chr(j)
      rcov2(i)=rcov2(i)+(chr(j)**2)
      end do
      bp=float(rend(i)-rbeg(i)+1)
      acov(i)=float(rcov(i))/bp
      acov2(i)=float(rcov2(i))-((float(rcov(i))/bp)*float(rcov(i)))
      acov2(i)=acov2(i)/(bp-1)
      end do
      open(unit=1,file='GC_'//arr//'_Coverage_Summary_Stats.txt')
      write(1,'("Gene Conversion Tract Summary Statistics ",a2)') arr
      write(1,'(1x)')
      write(1,'("Region  Mean_GC  Var__GC  SD___GC  X - 2SD  X + 2SD")') 
      do i=1,nr
      write(1,'(\i6,3f9.3)') i,acov(i),acov2(i),(acov2(i)**0.5)
      if(acov(i)-(2.*(acov2(i)**0.5)).lt.0.) then
      write(1,'(4x,"0.000",f9.3)') acov(i)+(2.*(acov2(i)**0.5))
      else
      write(1,'(\f9.3)') acov(i)-(2.*(acov2(i)**0.5))
      write(1,'(f9.3)') acov(i)+(2.*(acov2(i)**0.5))
      end if
      end do
      close(unit=1)
c
c-----Estimate Mean Gene Conversion Tract Coverage per CDS
c
c-----ggc  - Count of gene conversion tracts
c-----ggcs - Count of nucleotide sites
c
      gcov=0
      gsit=0.
      do i=1,ntr
      do j=ntb(i),nte(i)
      gcov(ts(i),rchr(j))=gcov(ts(i),rchr(j))+chr(j)
      gsit(ts(i),rchr(j))=gsit(ts(i),rchr(j))+1.
      end do
      end do
      do i=1,nr
      bcov(1,i)=float(gcov(1,i))/gsit(1,i)
      bcov(2,i)=float(gcov(2,i))/gsit(2,i)
      end do
      do i=1,nr      
      end do
      close(unit=1)
c      open(unit=1,file='GC_'//arr//'_Noncoding_List.tsv')
c      ntr=0
c350   ntr=ntr+1
c
c-----rname(i) - Name of the ith noncoding region
c-----ntex(i)  - Type of region [1=Intergenic; 2=Intron]
c-----ntb(i)   - First nucleotide of the ith transcript exon
c-----nte(i)   - Last nucleotide of the ith transcript exon
c-----ts(i)    - Status of the ith transcript exon [1=non-outlier; 2=outlier]
c
c
c-----GC_XX_Transcripts_Exons_Sort_List.tsv has a sorted list of exons
c-----for each gene.  The genes are listed in alphabetical order.
c
      open(unit=1,file='GC_'//arr//'_Transcripts_Exons_Sort_List.tsv')
      ntr=0
250   ntr=ntr+1
      read(1,*,end=300) gname(ntr),ntex(ntr),ntb(ntr),nte(ntr),ts(ntr)
      goto 250
300   close(unit=1)
      ntr=ntr-1
c
c-----Estimate GC Coverage in CDS
c-----Output to file GC_XX_CDS_Coverage.tsv
c
      open(unit=1,file='GC_'//arr//'_CDS_Coverage.tsv')
      write(1,'(\"Gene_Name",a1,"Beg",a1,"End",a1,"Region")') t,t,t
      write(1,'(a1,"Nuc",a1,"Mean_Cov",a1,"Outlier_Status")') t,t,t
      gnamt=gname(1)
      ggc=0
      ggcs=0
      imin=ntb(1)
      imax=nte(1)
      do i=1,ntr
      if(gname(i).eq.gnamt) then
      if(ntb(i).lt.imin) imin=ntb(i)
      if(nte(i).gt.imax) imax=nte(i)
      do j=ntb(i),nte(i)
      ggc=ggc+chr(j)
      ggcs=ggcs+1
      end do
      else
      write(1,'(\a11,a1,i10,a1,i10)') gname(i-1),t,imin,t,imax
      write(1,'(\a1,i2,a1,i10)') t,rchr(imin),t,ggcs
      write(1,'(a1,f9.3,a1,i1)') t,float(ggc)/float(ggcs),t,ts(i-1)
      gnamt=gname(i)
      ggc=0
      ggcs=0
      imin=ntb(i)
      imax=nte(i)
      do j=ntb(i),nte(i)
      ggc=ggc+chr(j)
      ggcs=ggcs+1
      end do
      end if
      end do
      write(1,'(\a11,a1,i10,a1,i10)') gname(ntr),t,imin,t,imax
      write(1,'(\a1,i2,a1,i10)') t,rchr(imin),t,ggcs
      write(1,'(a1,f9.3,a1,i1)') t,float(ggc)/float(ggcs),t,ts(ntr)
      close(unit=1)
c
c-----Extract intron coordinates
c
c-----Output to file GC_XX_Intron_List.tsv
c
      open(unit=1,file='GC_'//arr//'_Intron_List.tsv')
      write(1,'(\"Transcript",a1,"Intron_No",a1,"Beg",a1,"End")') t,t,t
      write(1,'(a1,"Outlier")')
      nin=0
      gnamt=gname(1)
      do i=2,ntr
      if(gnamt.eq.gname(i)) then
      nin=nin+1
      write(1,'(\a11,a1,i2)') gnamt,t,nin
      write(1,'(a1,i8,a1,i8,a1,i1)') t,nte(i-1)+1,t,ntb(i)-1,t,ts(i)
      else
      nin=0
      gnamt=gname(i)
      end if
      end do
      close(unit=1)
c
c-----Open intron list file and input the coordinates
c
      open(unit=1,file='GC_'//arr//'_Intron_List.tsv')
      read(1,'(1x)')
      ntr=0
350   ntr=ntr+1
      read(1,*,end=400) gname(ntr),ntex(ntr),ntb(ntr),nte(ntr),ts(ntr)
      goto 350
400   close(unit=1)
      ntr=ntr-1
c
c-----Estimate GC coverage in introns
c-----Output to GC_XX_Intron_Coverage.tsv
c      
      open(unit=1,file='GC_'//arr//'_Intron_Coverage.tsv')
      write(1,'(\"Gene_Name",a1,"Beg",a1,"End",a1,"Region")') t,t,t
      write(1,'(a1,"Nuc",a1,"Mean_Cov",a1,"Outlier_Status")') t,t,t
      gnamt=gname(1)
      ggc=0
      ggcs=0
      imin=ntb(1)
      imax=nte(1)
      do i=1,ntr
      if(gname(i).eq.gnamt) then
      if(ntb(i).lt.imin) imin=ntb(i)
      if(nte(i).gt.imax) imax=nte(i)
      do j=ntb(i),nte(i)
      ggc=ggc+chr(j)
      ggcs=ggcs+1
      end do
      else
      write(1,'(\a11,a1,i10,a1,i10)') gname(i-1),t,imin,t,imax
      write(1,'(\a1,i2,a1,i10)') t,rchr(imin),t,ggcs
      write(1,'(a1,f9.3,a1,i1)') t,float(ggc)/float(ggcs),t,ts(i-1)
      gnamt=gname(i)
      ggc=0
      ggcs=0
      imin=ntb(i)
      imax=nte(i)
      do j=ntb(i),nte(i)
      ggc=ggc+chr(j)
      ggcs=ggcs+1
      end do
      end if
      end do
      write(1,'(\a11,a1,i10,a1,i10)') gname(ntr),t,imin,t,imax
      write(1,'(\a1,i2,a1,i10)') t,rchr(imin),t,ggcs
      write(1,'(a1,f9.3,a1,i1)') t,float(ggc)/float(ggcs),t,ts(ntr) 
      close(unit=1)
c
c-----Use the CDS coordinates to extract the Intergenic Regions
c
c-----igchr(i) - array whose elements map CDS coordinates to it
c-----gamt     - name of the CDS transcript
c-----nib      - the beginning coordinate of the CDS     
c-----nie      - the end coordinate of the CDS
c
c-----Each CDS interval is used to increment the igchr array
c-----At the end, an element value of zero will be an intergenic nucleotide 
c
      igchr=0
      open(unit=1,file='GC_'//arr//'_Transcripts_CDS_List.tsv')
      read(1,'(1x)')
450   read(1,*,end=500) gnamt,nib,nie
c
c-----Map the CDS coordinates
c
      do i=nib,nie
      igchr(i)=igchr(i)+1
      end do
      goto 450 
500   close(unit=1)
      open(unit=1,file='GC_'//arr//'_Intergenic_List.tsv')
      write(1,'("IG_Region",a1,"Beg",a1,"End")') t,t
c
c-----Initiate variables
c-----nr    - number of the region
c-----nib   - beginning coordinate of the intergenic region
c-----nie   - end coordinate of the intergenic region
c-----iflag - variable used to keep track of intergenic or CDS intervals [0, intergenic; 1, CDS] 
c 
      nr=1
      nib=1
      nie=1
      iflag=0
      do i=1,19787792
      if(igchr(i).eq.0.and.iflag.eq.0) then
c
c-----Nucleotide is intergenic and already in an intergenic region iflag value (0)
c-----Increment nie 
c
      nie=i
      elseif(igchr(i).ne.0.and.iflag.eq.0) then
c
c-----Nucleotide is CDS and transition to CDS
c-----Switch iflag to CDS value (1)
c-----Output the intergenic interval
c
      iflag=1
      write(1,'("Region_",i4.4,a1,i8,a1,i8)') nr,t,nib,t,nie
      elseif(igchr(i).eq.0.and.iflag.eq.1) then
c
c-----Nucleotide is intergenic and transition to intergenic
c-----Switch iflag to intergenic value (0)
c-----Increment nr, Reset nib and nie
c
      iflag=0
      nr=nr+1
      nib=i
      nie=i 
      end if
      end do
c
c-----Output the last intergenic interval
c
      write(1,'("Region_",i4.4,a1,i8,a1,i8)') nr,t,nib,t,nie
      close(unit=1)
c
c-----Reopen file with intergenic coordinates and input regions and coordinates
c
c-----Reuse variables gname(i), ntb(i), and nte(i)
c-----Reuse variables ggc (coverage), ggcs (sites)
c
      open(unit=1,file='GC_'//arr//'_Intergenic_List.tsv')
      read(1,'(1x)')
      ntr=0
550   ntr=ntr+1
      read(1,*,end=600) gname(ntr),ntb(ntr),nte(ntr)
      goto 550
600   close(unit=1)
      ntr=ntr-1
c
c-----Estimate gene conversion coverage in intergenic regions
c
      open(unit=1,file='GC_'//arr//'_Intergenic_Coverage.tsv')
      write(1,'(\"Gene_Name",a1,"Beg",a1,"End",a1,"Region")') t,t,t
      write(1,'(a1,"Nuc",a1,"Mean_Cov")') t,t   
      do i=1,ntr
      ggc=0
      ggcs=0
      do j=ntb(i),nte(i)
      ggc=ggc+chr(j)
      ggcs=ggcs+1
      end do
c
c-----Output the results for the intergenic interval
c
      write(1,'(\a11,a1,i10,a1,i10)') gname(i),t,ntb(i),t,nte(i)
      write(1,'(\a1,i2,a1,i10)') t,rchr(ntb(i)),t,ggcs
      write(1,'(a1,f9.3)') t,float(ggc)/float(ggcs)
      end do      
      close(unit=1)
      stop
      end
