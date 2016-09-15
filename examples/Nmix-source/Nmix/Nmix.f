	parameter(nmax=1000,ncmax=30,ncmax2=10,ngrid=200,ncd=7,
     &		k1=2,k2=8,kkzz=((k2-k1+1)*(k1+k2))/2)

cd.. output file switches 
cd  standard: bdlog ent k log out pe
cd  extra with -debug: db
cd  extra with -beta (=default), -kappa: bk
cd  extra with -p*d* or -full: avn den dev
cd  extra with -p*c*: pcl scl
cd  extra with -p*w*: wms.*
cd  extra with -p*a*: z.*
cd
cd.. options
cd  Setting up run 
cd		default
cd  -n		10000
cd  -nb		0
cd  -ki		1
cd  -seed	random
cd  
cd  Model options
cd  (default: all off)
cd  -prior
cd  -fix
cd  -unif
cd  
cd  Parameter settings
cd		default
cd  -lu		on
cd  -l		0 
cd  -a		2
cd  -b		Gamma(g,h) (0.02 if fixed)
cd  -x 		0
cd  -k		1
cd  -unif		1
cd  -sp		1
cd  -d		1
cd  -e		0
cd  -f		0
cd  -g		0.2
cd  -h		10
cd  -range	on
cd  
cd  Output options
cd		default
cd  -full	off
cd  -p*d*	off
cd  -p*a*	off
cd  -p*w*	off
cd  -p*c*	off
cd  -nsp	20
cd  -ns	100
cd  -nsk	del2
cd  -debug	off
cd  -rkp	off
cd  
cd  Sampler switches
cd  (default: metro, empty on, others off)
cd  -metro 
cd  -empty	
cd  -o		5
cd  -m		swpahb for std

	character*18 word,base,prfx
	character*30 wordl
	character*6 num,num2
	character*80 line
	character moveselect
	character*10 moves
	character*5 movnam(8)
	character sep

	logical fileopen(ncmax)
	logical qbeta,qkappa,qq,qdebug,qprior,qrange,
     &		qfull,qfix,
     &		qunif,qempty,qstop,qpalloc,qpwms,qpclass,
     &		qkreg,qrkpos

	integer hwm,first,free,st1,st2,sp,ckrep
	integer prev(ncmax),next(ncmax),start(ncmax),leng(ncmax)
	integer inext(nmax),na(nmax),z(nmax)
	integer nfkemp(0:9),count(ncmax),countpos(ncmax),countden(ncd)
        integer iseed

	real kappa,lambda,mu1,mu2,muc,ms,logratio,logbeta,lprob
	real kwas,loglikrat,mustar
	real y(nmax),yv(nmax)
	real wt(ncmax),mu(ncmax),ssq(ncmax),lp(ncmax),mun(ncmax)
	real b(ncmax),d(ncmax)
	real wtav(ncmax2,ncmax2),muav(ncmax2,ncmax2)
	real sigav(ncmax2,ncmax2),avn(ncmax2,ncmax2)
	real den(ncd,ngrid),pclass(kkzz,ngrid),pz(kkzz,nmax)
	real pw(ncmax),bf(ncmax),avdev(ncmax),avdevc(ncmax),ggh(ncd)
	real pf(ncmax)
	real*8 dlgama
	integer split,combine,allocate,weights,parameters,hyper,
     &		birth,death
	data split,combine,allocate,weights,parameters,hyper,
     &		birth,death/1,2,3,4,5,6,7,8/
	data movnam/'split','comb ','alloc','wts  ','param','hyper',
     &		'birth','death'/
c.. subprograms:
c	sdrni(0)		random number initialiser
c	sdrand()		uniform random number
c	gauss(z,n)	z = standard normal random vector of length n
c	rgamma(alpha)	gamma random number, shape alpha, scale 1
c	rbeta(alpha1,alpha2) beta random number, parameters alpha1,alpha2

c.. notation for partitions:
c	Component j has weight wt(j), and parameters mu(j) and ssq(j); it
c       contains leng(j) observations y(i), i=i1=start(j),i2=inext(i1),i3=inext(i2),...
c	The component with the next smallest mu is prev(j), and
c       that with next largest mu is next(j); in either case if there
c       is no such mu, the pointer is set to 0.
c       Non-existent components have prev(j) = -1.

	sep = "/"
	pi = 4.0*atan(1.0)
	rr2pi = 1.0/sqrt(2.0*pi)

c.. model hyperparameters
c	each mu is N(xi,kappa^{-1})
c	each 1/ssq is Gamma(alpha,beta)
c	weights are Dirichlet(delta,...,delta)
c	k is Poisson(lambda), conditioned to be positive
c	  lambda = 0 => p(k) = 1/k
c       lambda = -1 => p(k) = uniform
c	beta is gamma(g,h), or fixed
c	kappa is gamma(e,f), or fixed

	qstop = .false.

	xi = 0.0
	kappa = 1.0
	qkappa = .false.
	unhw = 1.0
	eee = 0.0
	fff = 0.0
	sp = 1
	alpha = 2.0
	beta = 0.02
	qbeta = .true.
	ggg = 0.2
	hhh = 10.0
	delta = 1.0
	lambda = -1.0
	omega = 5.0
	kinit = 1

	nsweep = 10000
	nburnin = 0
	nsamp = 100
	nspace = 20
	nskdel = 2
	idebug = -1
	qdebug = .false.
	qunif = .false.
	qempty = .true.
	qprior = .false.
	qrange = .true.
	qfull = .false.
	qpalloc = .false.
	qpwms = .false.
	qpclass = .false.
	qkreg = .false.
	qrkpos = .false.
	qfix = .false.

	iseed = 0
c.. pick up base for filenames, removing ".dat" if necessary

	if(iargc().lt.1) then
		stop 77
	end if
	do ia = 1,iargc()
		call getarg(ia,word)
		if(word(1:1).ne.'-') go to 2
	end do
	stop 76

2	base = word

	do ncbase = 18,1,-1
		if(base(ncbase:ncbase).ne." ") go to 7
	end do
	stop 75

7	if(ncbase.ge.5) then
	    if(base(ncbase-3:ncbase).eq.".dat") then
		base = base(1:ncbase-4)
		ncbase = ncbase-4
	    end if
	end if

	inquire(file=base(1:ncbase)//".dat",exist=qq)
	if(.not.qq) then
		write(0,'("data file not there:",a18)')
     &		base(1:ncbase)//".dat"
		stop
	end if

c.. read control parameters, from files first if present
c.. copy command name for log file

	iline = 0
	call getarg(0,wordl)
	do j = 1,30
	    iline = iline+1
	    line(iline:iline) = wordl(j:j)
	    if(wordl(j:j).eq.' ') go to 54
	end do
54	continue

	ia = -3
	do while(ia.le.iargc())
	if(ia.eq.-3) then
		inquire(file="pars",exist=qq)
		if(qq) then
			open(3,file="pars",status='unknown')
			ia = -2
		else
			ia = -1
		end if
		go to 50
	else if(ia.eq.-2) then
		read(3,'(a18)',end=51) word
		go to 52
51		close(3)
		ia = -1
		go to 50
	else if(ia.eq.-1) then
		inquire(file=base(1:ncbase)//".pars",exist=qq)

		if(qq) then
			open(3,file=base(1:ncbase)//".pars",status='unknown')

			ia = 0
		else
			ia = 1
		end if
		go to 50
	else if(ia.eq.0) then
		read(3,'(a18)',end=53) word
		go to 52
53		close(3)
		ia = 1
		go to 50	
	else
		call getarg(ia,word)
		ia = ia+1
	end if
52	if(word(1:3).eq.'-nb') then
		read(word,'(3x,i15)') nburnin
	else if(word(1:4).eq.'-nsp') then
		read(word,'(4x,i14)') nspace
	else if(word(1:4).eq.'-nsk') then
		if(word(5:7).eq.'del') then
			read(word,'(7x,i11)') nskdel 
		else
			qkreg = .true.
			read(word,'(4x,i14)') nskdel
		end if
	else if(word(1:3).eq.'-ns') then
		read(word,'(3x,i15)') nsamp
	else if(word(1:3).eq.'-ki') then
		read(word,'(3x,i15)') kinit
		kinit = max(0,min(ncmax,kinit))
	else if(word(1:5).eq.'-seed') then
		read(word,'(5x,i13)') iseed
	else if(word(1:4).eq.'-fix') then
		qfix = .true.
	else if(word(1:5).eq.'-unif') then
		qunif = .true.
		read(word,'(5x,f13.0)') unhwin
		if(unhwin.gt.0.0) unhw = unhwin
	else if(word(1:6).eq.'-range') then
		qrange = .true.
	else if(word(1:8).eq.'-norange') then
		qrange = .false.
	else if(word(1:6).eq.'-empty') then
		qempty = .true.
	else if(word(1:8).eq.'-noempty') then
		qempty = .false.
	else if(word(1:2).eq.'-n') then
		read(word,'(2x,i16)') nsweep
	else if(word(1:6).eq.'-debug') then
		read(word,'(6x,i12)') idebug
	else if(word(1:6).eq.'-prior') then
		qprior = .true.
	else if(word(1:5).eq.'-full') then
		qfull = .true.
	else if(word(1:2).eq.'-p') then
		do j = 3,18
		    if(word(j:j).eq.'d') then
			qfull = .true.
		    else if(word(j:j).eq.'a') then
			qfull = .true.
			qpalloc = .true.
		    else if(word(j:j).eq.'w') then
			qfull = .true.
			qpwms = .true.
		    else if(word(j:j).eq.'c') then
			qfull = .true.
			qpclass = .true.
		    end if
		end do
	else if(word(1:4).eq.'-rkp') then
		qrkpos = .true.
	else if(word(1:3).eq.'-lu') then
		lambda = -1.0
	else if(word(1:2).eq.'-l') then
		read(word,'(2x,f16.0)') lambda
	else if(word(1:2).eq.'-a') then
		read(word,'(2x,f16.0)') alpha
	else if(word(1:2).eq.'-b') then
		read(word,'(2x,f16.0)') beta
		if(beta.le.0.0) then
			qbeta = .true.
			beta = 0.02
		else
			qbeta = .false.
		end if
	else if(word(1:2).eq.'-x') then
		read(word,'(2x,f16.0)') xi
	else if(word(1:2).eq.'-k') then
		read(word,'(2x,f16.0)') kappa
		if(kappa.le.0.0) then
			qkappa = .true.
			kappa = 4.0
		else
			qkappa = .false.
		end if
	else if(word(1:3).eq.'-sp') then
		read(word,'(3x,i15)') sp
		sp = max(1,sp)
	else if(word(1:2).eq.'-o') then
		read(word,'(2x,f16.0)') omega
	else if(word(1:2).eq.'-d') then
		read(word,'(2x,f16.0)') delta
	else if(word(1:2).eq.'-e') then
		read(word,'(2x,f16.0)') eee
	else if(word(1:2).eq.'-f') then
		read(word,'(2x,f16.0)') fff
	else if(word(1:2).eq.'-g') then
		read(word,'(2x,f16.0)') ggg
	else if(word(1:2).eq.'-h') then
		read(word,'(2x,f16.0)') hhh
	end if

c.. copy argument for log file

	do j = 1,18
	    iline = iline+1
	    line(iline:iline) = word(j:j)
	    if(word(j:j).eq.' ') go to 55
	end do
55	continue

50	continue

	end do

c.. final parameter fixing

	call sdrni(iseed)

	if(qunif) then
		sp = 1
		qkappa = .false.
	end if

	qempty = qempty.and.(sp.eq.1).and.(.not.qfix)

	if(nburnin.lt.0) nburnin = nsweep/10
	if(idebug.lt.0) idebug = nsweep+1
	qdebug = idebug.eq.0

	if(qdebug) then
		nburnin = 0
		nsamp = nsweep
	end if

c.. sweep schedule

c standard: sequence is 0 1 2 3 4
c standard/fix: sequence is 1 2 3 4


	if(qfix) then
	   moves = 'wpah'
	   nstep = 4
	else
	   moves = 'swpah'
	   nstep = 5
	end if

	if(qempty) then
		nstep = nstep+1
		moves(nstep:nstep) = 'b'
	end if


c.. sort out filename stuff

	ijk = 1

1	istd = 6-int(log10(0.5+ijk))
	write(num,'(i6)') ijk

	inquire(file=base(1:ncbase)//sep//num(istd:6)//".log",exist=qq)
	if(qq) then
		ijk = ijk+1
		go to 1
	end if

	prfx = base(1:ncbase)//sep//num(istd:6)
	npf = ncbase+8-istd

	open(6,file=prfx(1:npf)//".out",status='unknown')
	open(7,file=prfx(1:npf)//".log",status='unknown')
	if(.not.qfix) open(8,file=prfx(1:npf)//".bdlog",status='unknown')

	if(qfull) then
	    open(9,file=prfx(1:npf)//".dev",status='unknown')
	    if(qpclass) then
		open(14,file=prfx(1:npf)//".pcl",status='unknown')
		open(15,file=prfx(1:npf)//".scl",status='unknown')
	    end if
	end if
	if(qbeta.or.qkappa)
     &	    open(10,file=prfx(1:npf)//".bk",status='unknown')
	open(12,file=prfx(1:npf)//".k",status='unknown')
	if(qfull)
     &		open(13,file=prfx(1:npf)//".den",status='unknown')
	open(17,file=prfx(1:npf)//".ent",status='unknown')
	open(11,file=prfx(1:npf)//".pe",status='unknown')
	if(qfull) then
	    open(16,file=prfx(1:npf)//".avn",status='unknown')
	end if

	if(idebug.le.nsweep) then
	    open(18,file=prfx(1:npf)//".db",status='unknown')
	end if

c.. proposal parameters

	ws = 2.0
	ms = 2.0
	ss = 1.0

	do k = 2,ncmax-1
		b(k) = 0.5
		d(k) = 0.5
	end do
	b(ncmax) = 0.0
	d(ncmax) = 1.0
	b(1) = 1.0
	d(1) = 0.0

c.. read data, and compute basic summary statistics

	open(5,file=base(1:ncbase)//".dat",status='unknown')
	read(5,*) n
	if(n.gt.nmax) stop 98
	read(5,*) (y(i),i=1,n)
	close(5)

c.. basic statistics

	ymin = y(1)
	ymax = ymin
	do i = 2,n
		ymin = min(ymin,y(i))
		ymax = max(ymax,y(i))
	end do

	ymid = 0.5*(ymin+ymax)
	yrange = ymax-ymin
	yd0 = ymin-0.05*yrange
	ydinc = 1.1*yrange/ngrid

	ysum = 0.0
	do i = 1,n
		ysum = ysum+y(i)
	end do
	ysum = ysum/n
	ssd  = 0.0
	do i = 1,n
		ssd = ssd+(y(i)-ysum)**2
	end do
c	devcorrection = -n*log(ssd/(n-1))

c.. adjust hyperparameters if range-based

	if(qrange) then
		xiwas = xi
		xi = ymid+xi*0.5*yrange
		kwas = kappa
		kappa = kappa/yrange**2
		unhwwas = unhw
		unhw = unhw*0.5*yrange
		bwas = beta
		beta = beta*yrange**2
		fwas = fff
		fff = fff*yrange**2
		hwas = hhh
		hhh = hhh/yrange**2
	end if

c.. set up prior on k

	    if(lambda.gt.0.0) then
		temp = -log(exp(lambda)-1.0)
		do k = 1,ncmax
			temp = temp+log(lambda/k)
			lp(k) = temp
		end do
	    else if(lambda.lt.0.0) then
		do k = 1,ncmax
			lp(k) = 0.0
		end do
	    else
		do k = 1,ncmax
			lp(k) = log(1.0/k)
		end do
	    end if

c.. log parameter values

	write(0,*) "Run: ",prfx(1:npf)," : ",line(1:iline)
	write(7,*) "Run: ",prfx(1:npf)," : ",line(1:iline)
	write(7,'("random number seed:",i12)') iseed
	if(qprior) write(7,'("PRIOR simulation")')
	if(qunif) write(7,'("Uniform prior for component means")')
	if(qfix) write(7,'("fixed k")')
	if(qfull) write(7,'("full output")')
	if(qrkpos) write(7,'("reporting kpos option")')
	if(qrange) then
		write(7,
     &	 '("xi      =   ",g11.4,"*R/2+ymid =",g11.4)') xiwas,xi
	else
		write(7,'("xi      =   ",g11.4)') xi
	end if
	if(qunif) then
		if(qrange) then
			write(7,
     &	'(" unhw = ",g11.4,"*R/2 =",g11.4)') unhwwas,unhw
		else
			write(7,'(" unhw = ",g11.4)') unhw
		end if
	else if(qkappa) then
		write(7,'("kappa is Gamma(e,f) with")')
		write(7,'("  e     =   ",g11.4)') eee
		if(qrange) then
			write(7,
     &	 '("  f     =   ",g11.4,"*R^2 =",g11.4)') fwas,fff
		else
			write(7,'("  f     =   ",g11.4)') fff
		end if
	else
		if(qrange) then
			write(7,
     &	 '("kappa   =   ",g11.4,"/R^2 =",g11.4)') kwas,kappa
		else
			write(7,'("kappa   =   ",g11.4)') kappa
		end if
	end if
	if(sp.gt.1) write(7,'(" spacing prior, sp=",i3)') sp
	write(7,'("alpha   =   ",g11.4)') alpha
	if(qbeta) then
		write(7,'("beta is Gamma(g,h) with")')
		write(7,'("  g     =   ",g11.4)') ggg
		if(qrange) then
			write(7,
     &	 '("  h     =   ",g11.4,"/R^2 =",g11.4)') hwas,hhh
		else
			write(7,'("  h     =   ",g11.4)') hhh
		end if
	else
		if(qrange) then
			write(7,
     &	 '("beta   =   ",g11.4,"*R^2 =",g11.4)') bwas,beta
		else
			write(7,'("beta    =   ",g11.4)') beta
		end if
	end if
	write(7,'("delta   =   ",g11.4)') delta
	write(7,'("lambda  =   ",g11.4)') lambda
	write(7,'("nsweep  =   ",i8)') nsweep
	write(7,'("nstep   =   ",i8)') nstep
	write(7,'("nburnin =   ",i8)') nburnin
	write(7,'("n       =   ",i8)') n
	write(7,'("ymin    =   ",g11.4)') ymin
	write(7,'("ymax    =   ",g11.4)') ymax

	write(7,*) "move schedule: ",moves(1:nstep)

c.. initialise dynamic variables
c.. initialise number of components, k, and allocation

	if(kinit.eq.0) then
		open(3,file=base(1:ncbase)//".zin",status='unknown')
		read(3,*) (na(i),i=1,n)
		close(3)
		do i = 1,n
			if(na(i).lt.1.or.na(i).gt.ncmax) stop 74
		end do
		do j = 1,ncmax
			start(j) = 0
			leng(j) = 0
		end do
		do i = 1,n
			j = na(i)
			inext(i) = start(j)
			start(j) = i
			leng(j) = leng(j)+1
		end do
		k = 0
		first = 0
		free = 0
		do j = 1,ncmax
		    next(j) = 0
		    if(leng(j).gt.0) then
			hwm = j
			k = k+1
			if(first.eq.0) then
				first = j
				prev(j) = 0
			else
				prev(j) = jp
				next(jp) = j
			end if
			jp = j
		    else
			if(free.eq.0) then
				free = j
			else
				next(jq) = j
			end if
			prev(j) = -1
			jq = j
		    end if
		end do
	write(7,'("initial allocation read from .zin file: k =",i3)') k
	else
	write(7,'("initial number of components:",i4)') kinit
	if(kinit.eq.1) then
	    start(1) = 1
	    leng(1) = n
	    do i = 1,n
		inext(i) = i+1
	    end do
	    inext(n) = 0
	else
	    do i = 1,n
		na(i) = i
	    end do
	    call sort(na,1,n,y,n)
	    do j = 1,kinit
		l1 = (n*(j-1))/kinit+1
		l2 = (n*j)/kinit
		leng(j) = l2-l1+1
		if(leng(j).gt.0) then
			start(j) = na(l1)
		else
			start(j) = 0
		end if
		do l = l1,l2-1
			inext(na(l)) = na(l+1)
		end do
		if(l2.gt.0) inext(na(l2)) = 0
	    end do
	end if
	    first = 1
	    do j = 1,kinit
		prev(j) = j-1
		next(j) = j+1
	    end do
	    next(kinit) = 0
	    free = kinit+1
	    do j = kinit+1,ncmax
		prev(j) = -1
		next(j) = j+1
	    end do
	    next(ncmax) = 0
	    k = kinit
	    hwm = kinit
	end if

	kemp = 0
	j = first
	do while(j.ne.0)
	    if(leng(j).eq.0) kemp = kemp+1
	    j = next(j)
	end do

c.. initialise weights

	wtsum = 0.0
	j = first
	do while(j.ne.0)
	    wt(j) = rgamma(delta+leng(j))
	    wtsum = wtsum+wt(j)
	    j = next(j)
	end do
	j = first
	do while(j.ne.0)
	    wt(j) = wt(j)/wtsum
	    j = next(j)
	end do

c.. initialise means

	call gauss(mu,hwm)
	if(qdebug) then
		write(18,*) (mu(j),j=1,hwm)
	end if
	ssq0 = beta/alpha
	j = first
	do while(j.ne.0)
	if(qprior) then
	    if(qunif) then
		mu(j) = xi+unhw*(2.0*sdrand()-1.0)
	    else
		mu(j) = xi+mu(j)/sqrt(kappa)
	    end if
	else
	    ysum = 0.0
	    i = start(j)
	    do while(i.ne.0)
		ysum = ysum+y(i)
		i = inext(i)
	    end do
	    if(qunif) then
		if(leng(j).eq.0) then
		    mu(j) = xi+unhw*(2.0*sdrand()-1.0)
		else
44		    mu(j) = ysum/leng(j)+sqrt(ssq0/leng(j))*mu(j)
		    if(abs((mu(j)-xi)/unhw).gt.1.0) then
			call gauss(mu(j),1)
			go to 44
		    end if
		end if
	    else
		con = 1.0/(leng(j)/ssq0+kappa)
		mu(j) = con*(ysum/ssq0+xi*kappa)+sqrt(con)*mu(j)
	    end if
	end if
	    j = next(j)
	end do

c.. initialise variances

	j = first
	do while(j.ne.0)
	if(qprior) then
	    ssq(j) = beta/rgamma(alpha)
	else
	    ssd = 0.0
	    i = start(j)
	    do while(i.ne.0)
		ssd = ssd+(y(i)-mu(j))**2
		i = inext(i)
	    end do
	    ssq(j) = (beta+0.5*ssd)/rgamma(alpha+0.5*leng(j))
	end if
	j = next(j)
	end do

c.. now check ordering of mu, and correct

	call reorder(mu,ncmax,next,prev,first)

c.. initialise v
	do i = 1,n
		yv(i) = y(i)
	end do

c.. initialise accumulators

	do j = 1,ncmax
		count(j) = 0
		countpos(j) = 0
		fileopen(j) = .false.
		avdev(j) = 0.0
		avdevc(j) = 0.0
	end do

	do j = 1,ncmax2
	do ij = 1,j
		wtav(j,ij) = 0.0
		muav(j,ij) = 0.0
		sigav(j,ij) = 0.0
		avn(j,ij) = 0.0
	end do
	end do

	do krec = 1,ncd
	    countden(krec) = 0
	    do iyd = 1,ngrid
		den(krec,iyd) = 0.0
	    end do
	end do

	nkemp = 0
	avkemp = 0.0
	ppkemp = 0.0
	do j = 0,9
		nfkemp(j) = 0
	end do

	if(qpclass) then
	    do iyd = 1,ngrid
	    do j = 1,kkzz
		pclass(j,iyd) = 0.0
	    end do
	    end do

	    do i = 1,n
	    do j = 1,kkzz
		pz(j,i) = 0.0
	    end do
	    end do
	end if

	ntrys = 0
	ntryc = 0
	ntryb = 0
	ntryd = 0
	naccs = 0
	naccc = 0
	naccb = 0
	naccd = 0
	nrejr = 0

	if(qunif) then
	    const1 = -sngl(dlgama(dble(alpha)))
	else
	    const1 = 0.5*log(0.5/pi)-sngl(dlgama(dble(alpha)))
	end if
	const2 = logbeta(ws,ws)+logbeta(ms,ms)+logbeta(ss,ss)
	if(qdebug) write(18,'("consts",2g11.5)') const1,const2

	do j = 1,ncmax
	temp = sp*(j+1)
	    do jj = 1,sp-1
		temp = temp*real(sp*(j+1)+jj)/real(jj)
	    end do
	pf(j) = temp
	end do

	if(qdebug) then
		write(18,'("initialised")')
		write(18,*) k
		j = first
		do while(j.ne.0)
		    write(18,'(i4,i6,3f10.4)') j,leng(j),mu(j),ssq(j),wt(j)
		    j = next(j)
		end do
	end if

c.. main loop: MCMC iteration

	do isweep = 1-nburnin,nsweep

	if(mod(nsweep-isweep,nsweep/10).eq.0) then
	    write(0,'(i3,$)') (nsweep-isweep)/(nsweep/10)
	end if

	qdebug = isweep.ge.idebug

	do istep = 0,nstep-1
c.. select next move type

	moveselect = moves(istep+1:istep+1)
	if(moveselect.eq.'s') then
	    usw = sdrand()
	    if(qdebug) write(18,*) isweep,usw
	    if(usw.lt.b(k)) then
		move = split
	    else
		move = combine
	    end if
	else if(moveselect.eq.'w') then
		move = weights
	else if(moveselect.eq.'p') then
		move = parameters
	else if(moveselect.eq.'a') then
		move = allocate
	else if(moveselect.eq.'h') then
		move = hyper
	else if(moveselect.eq.'b') then
	    usw = sdrand()
	    if(usw.lt.b(k)) then
		move = birth
	    else
		move = death
	    end if
	end if

	if(move.eq.split) then

c	   write(0,*)'split'

c.. split move ----------------------------------------------

	ntrys = ntrys+1

20	j = 1+int(hwm*sdrand())
	if(prev(j).lt.0) go to 20

	if(qdebug) write(18,*) isweep,j

	wtc = wt(j)
	muc = mu(j)
	ssqc = ssq(j)

	u1 = rbeta(ws,ws)
	cu1 = 1.0-u1
	u2 = rbeta(ms,ms)
	if(qdebug) write(18,*) isweep,u1,u2
	wt1 = wtc*u1
	wt2 = wtc-wt1
	mu1 = muc-sqrt(ssqc*wt1*wt2)*u2/wt1
	mu2 = (wtc*muc-wt1*mu1)/wt2

c.. check mu in order: else reject

	j1 = prev(j)
	if(j1.ne.0) then
		if(qdebug) write(18,*) mu1,mu(j1)
		if(mu(j1).gt.mu1) then
			nrejr = nrejr+1
			go to 24
		end if
	else if(qunif) then
		if(mu1.lt.xi-unhw) go to 24
	end if
	j1 = next(j)
	if(j1.ne.0) then
		if(qdebug) write(18,*) mu2,mu(j1)
		if(mu(j1).lt.mu2) then
			nrejr = nrejr+1
			go to 24
		end if
	else if(qunif) then
		if(mu2.gt.xi+unhw) go to 24
	end if

	u3 = rbeta(ss,ss)
	if(qdebug) write(18,*) isweep,u3
	cu3 = 1.0-u3
	temp = wtc*ssqc*(1.0-u2**2)
	ssq1 = temp*u3/wt1
	ssq2 = temp*cu3/wt2

c.. allocate observations in component to be split

	lprob = 0.0
	st1 = 0
	st2 = 0
	l1 = 0
	l2 = 0
	i = start(j)
	do while(i.ne.0)
	if(qprior) then
	    p1 = wt1
	    p2 = wt2
	else
	    p1 = wt1*exp(max(-20.0,-0.5*(yv(i)-mu1)**2/ssq1))/sqrt(ssq1)
	    p2 = wt2*exp(max(-20.0,-0.5*(yv(i)-mu2)**2/ssq2))/sqrt(ssq2)
	end if
	    in = inext(i)
	    if(sdrand().lt.p1/(p1+p2)) then
		inext(i) = st1
		if(st1.eq.0) ilast = i
		st1 = i
		l1 = l1+1
		lprob = lprob+log(p1)-log(p1+p2)
	    else
		inext(i) = st2
		st2 = i
		l2 = l2+1
		lprob = lprob+log(p2)-log(p1+p2)
	    end if
	    i = in
	end do

	if(qdebug) write(18,*) isweep,l1,l2

	klow = k

c.. compute ratio and decide

	logratio = 0.0
	if(.not.qprior) then
	i = st1
	do while(i.ne.0)
	    logratio = logratio-0.5*log(ssq1/ssqc)
     &		-0.5*((yv(i)-mu1)**2/ssq1-(yv(i)-muc)**2/ssqc)
	    i = inext(i)
	end do
	i = st2
	do while(i.ne.0)
	    logratio = logratio-0.5*log(ssq2/ssqc)
     &		-0.5*((yv(i)-mu2)**2/ssq2-(yv(i)-muc)**2/ssqc)
	    i = inext(i)
	end do
	end if

	if(qdebug) then
		write(18,'(2i6,1pe16.7)') isweep,move,logratio
	end if

c.. p(k,w,z) terms

	logratio = logratio+(lp(klow+1)-lp(klow))
     &		+(delta-1.0)*log(wt1*wt2/wtc)-logbeta(delta,klow*delta)
     &		+l1*log(wt1/wtc)+l2*log(wt2/wtc)

	if(qdebug) then
		write(18,'(2i6,1pe16.7)') isweep,move,logratio
	end if

c.. p(theta) terms

	if(qunif) then
	    logratio = logratio+const1+alpha*log(beta)
     &		-log(2.0*unhw)
     &		-(alpha+1.0)*log(ssq1*ssq2/ssqc)
     &		-beta*(1.0/ssq1+1.0/ssq2-1.0/ssqc)
     &		+log(pf(klow))
	else
	    logratio = logratio+const1+alpha*log(beta)+0.5*log(kappa)
     &		-0.5*kappa*((mu1-xi)**2+(mu2-xi)**2-(muc-xi)**2)
     &		-(alpha+1.0)*log(ssq1*ssq2/ssqc)
     &		-beta*(1.0/ssq1+1.0/ssq2-1.0/ssqc)
     &		+log(pf(klow))
	    if(sp.gt.1) then
		sqrk = sqrt(kappa)
		jplus = next(j)
		jminus = prev(j)
		if(jminus.eq.0) then
		    temp = log(pnorm((mu1-xi)*sqrk)/
     &			pnorm((muc-xi)*sqrk))
		else
		    temp = log((pnorm((mu1-xi)*sqrk)
     &			-pnorm((mu(jminus)-xi)*sqrk))
     &			/(pnorm((muc-xi)*sqrk)
     &			-pnorm((mu(jminus)-xi)*sqrk)))
		end if
		temp = temp+log(pnorm((mu2-xi)*sqrk)-pnorm((mu1-xi)*sqrk))
		if(jplus.eq.0) then
		    temp = temp+log((1.0-pnorm((mu2-xi)*sqrk))/
     &			(1.0-pnorm((muc-xi)*sqrk)))
		else
		    temp = temp+log((pnorm((mu(jplus)-xi)*sqrk)
     &			-pnorm((mu2-xi)*sqrk))
     &			/(pnorm((mu(jplus)-xi)*sqrk)
     &			-pnorm((muc-xi)*sqrk)))
		end if

		logratio = logratio+(sp-1)*temp
	    end if
	end if

	if(qdebug) then
		write(18,'(2i6,1pe16.7)') isweep,move,logratio
	end if

c.. proposal terms


	logratio = logratio+const2
     &		+log(d(klow+1)/b(klow))-lprob
     &		-(ws-1.0)*log(u1*cu1)
     &		-(ms-1.0)*log(u2*(1.0-u2))
     &		-(ss-1.0)*log(u3*cu3)

	if(qdebug) then
		write(18,'(2i6,1pe16.7)') isweep,move,logratio
	end if

c.. Jacobian terms

	logratio = logratio
     &		+log(wtc*abs(mu1-mu2)*ssq1*ssq2/ssqc)
     &		-log(u2*(1.0-u2**2)*u3*cu3)

	if(qdebug) then
		write(18,'(2i6,1pe16.7)') isweep,move,logratio
	end if

	logratio = max(-20.0,min(20.0,logratio))

      if(sdrand().lt.exp(logratio)) then

c.. accept split

		if(qdebug) then
		    write(18,'("split accepted")')
		end if
		naccs = naccs+1

		if(free.eq.0) stop 97
		jnew = free
		free = next(jnew)

		jnext = next(j)
		next(jnew) = jnext
		if(jnext.ne.0) prev(jnext) = jnew
		next(j) = jnew
		prev(jnew) = j

		hwm = max(hwm,jnew)

		start(j) = st1
		start(jnew) = st2
		leng(j) = l1
		leng(jnew) = l2
		wt(j) = wt1
		wt(jnew) = wt2
		mu(j) = mu1
		mu(jnew) = mu2
		ssq(j) = ssq1
		ssq(jnew) = ssq2

		k = k+1
		if(l1*l2.eq.0) kemp = kemp+1
		kpos = k-kemp

		if(.not.qkreg.and..not.qfix) write(8,'(2i3,i8)') k,kpos,isweep

	else

	    if(l1.ne.0) then
		start(j) = st1
		inext(ilast) = st2
	    else
		start(j) = st2
	    end if

	end if

24	continue

	else if(move.eq.combine) then

c	   write(0,*)'combine'

c.. combine move ----------------------------------------------

	ntryc = ntryc+1


30	j1 = 1+int(hwm*sdrand())
	if(prev(j1).lt.0) go to 30
	j2 = next(j1)
	if(j2.eq.0) go to 30


	st1 = start(j1)
	st2 = start(j2)
	l1 = leng(j1)
	l2 = leng(j2)

	wt1 = wt(j1)
	wt2 = wt(j2)
	mu1 = mu(j1)
	mu2 = mu(j2)
	ssq1 = ssq(j1)
	ssq2 = ssq(j2)
	wtc = wt1+wt2
	muc = (wt1*mu1+wt2*mu2)/wtc
c	ssqc = (wt1*(mu1**2+ssq1)
c     &		+wt2*(mu2**2+ssq2))/wtc-muc**2
	ssqc = wt1*wt2*((mu1-mu2)/wtc)**2 +(wt1*ssq1+wt2*ssq2)/wtc
	u1 = wt1/wtc
	cu1 = wt2/wtc
c	u2 = (muc-mu1)*wt1/sqrt(ssqc*wt1*wt2)
	u2 = (mu2-mu1)*sqrt(wt1*wt2/ssqc)/wtc
	u2 = max(u2,1e-12)
	u2 = min(u2,1.0-1e-4)
	u3 = wt1*ssq1/(wt1*ssq1+wt2*ssq2)
	cu3 = wt2*ssq2/(wt1*ssq1+wt2*ssq2)

c	write(0,'(3f8.5)')wt1,wt2,wtc,mu1,mu2,muc,ssq1,ssq2,ssqc,u1,u2,u3
c	write(0,*)ssq1,ssq2,ssqc,u2

	lprob = 0.0
	i = st1
	do while(i.ne.0)
	if(qprior) then
	    p1 = wt1
	    p2 = wt2
	else
	    p1 = wt1*exp(max(-20.0,-0.5*(yv(i)-mu1)**2/ssq1))/sqrt(ssq1)
	    p2 = wt2*exp(max(-20.0,-0.5*(yv(i)-mu2)**2/ssq2))/sqrt(ssq2)
	end if
	    lprob = lprob+log(p1)-log(p1+p2)
	    ilast = i
	    i = inext(i)
	end do
	i = st2
	do while(i.ne.0)
	if(qprior) then
	    p1 = wt1
	    p2 = wt2
	else
	    p1 = wt1*exp(max(-20.0,-0.5*(yv(i)-mu1)**2/ssq1))/sqrt(ssq1)
	    p2 = wt2*exp(max(-20.0,-0.5*(yv(i)-mu2)**2/ssq2))/sqrt(ssq2)
	end if
	    lprob = lprob+log(p2)-log(p1+p2)
	    i = inext(i)
	end do

	klow = k-1

c.. compute ratio and decide

	logratio = 0.0
	if(.not.qprior) then
	i = st1
	do while(i.ne.0)
	    logratio = logratio-0.5*log(ssq1/ssqc)
     &		-0.5*((yv(i)-mu1)**2/ssq1-(yv(i)-muc)**2/ssqc)
	    i = inext(i)
	end do
	i = st2
	do while(i.ne.0)
	    logratio = logratio-0.5*log(ssq2/ssqc)
     &		-0.5*((yv(i)-mu2)**2/ssq2-(yv(i)-muc)**2/ssqc)
	    i = inext(i)
	end do
	end if

	loglikrat = logratio

	if(qdebug) then
		write(18,'(2i6,1pe16.7)') isweep,move,logratio
	end if

c.. p(k,w,z) terms

	logratio = logratio+(lp(klow+1)-lp(klow))
     &		+(delta-1.0)*log(wt1*wt2/wtc)-logbeta(delta,klow*delta)
     &		+l1*log(wt1/wtc)+l2*log(wt2/wtc)

	if(qdebug) then
		write(18,'(2i6,1pe16.7)') isweep,move,logratio
	end if

c.. p(theta) terms

	if(qunif) then
	    logratio = logratio+const1+alpha*log(beta)
     &		-log(2.0*unhw)
     &		-(alpha+1.0)*log(ssq1*ssq2/ssqc)
     &		-beta*(1.0/ssq1+1.0/ssq2-1.0/ssqc)
     &		+log(pf(klow))
	else
	    logratio = logratio+const1+alpha*log(beta)+0.5*log(kappa)
     &		-0.5*kappa*((mu1-xi)**2+(mu2-xi)**2-(muc-xi)**2)
     &		-(alpha+1.0)*log(ssq1*ssq2/ssqc)
     &		-beta*(1.0/ssq1+1.0/ssq2-1.0/ssqc)
     &		+log(pf(klow))
	    if(sp.gt.1) then
		sqrk = sqrt(kappa)
		jplus = next(j2)
		jminus = prev(j1)
		if(jminus.eq.0) then
		    temp = log(pnorm((mu1-xi)*sqrk)/
     &			pnorm((muc-xi)*sqrk))
		else
		    temp = log((pnorm((mu1-xi)*sqrk)
     &			-pnorm((mu(jminus)-xi)*sqrk))
     &			/(pnorm((muc-xi)*sqrk)
     &			-pnorm((mu(jminus)-xi)*sqrk)))
		end if
		temp = temp+log(pnorm((mu2-xi)*sqrk)-pnorm((mu1-xi)*sqrk))
		if(jplus.eq.0) then
		    temp = temp+log((1.0-pnorm((mu2-xi)*sqrk))/
     &			(1.0-pnorm((muc-xi)*sqrk)))
		else
		    temp = temp+log((pnorm((mu(jplus)-xi)*sqrk)
     &			-pnorm((mu2-xi)*sqrk))
     &			/(pnorm((mu(jplus)-xi)*sqrk)
     &			-pnorm((muc-xi)*sqrk)))
		end if

		logratio = logratio+(sp-1)*temp
	    end if
	end if

	if(qdebug) then
		write(18,'(2i6,1pe16.7)') isweep,move,logratio
	end if


c.. proposal terms

	logratio = logratio+const2
     &		+log(d(klow+1)/b(klow))-lprob
     &		-(ws-1.0)*log(u1*cu1)
     &		-(ms-1.0)*log(u2*(1.0-u2))
     &		-(ss-1.0)*log(u3*cu3)

	if(qdebug) then
		write(18,'(2i6,1pe16.7)') isweep,move,logratio
	end if


c.. Jacobian terms

	logratio = logratio
     &		+log(wtc*abs(mu1-mu2)*ssq1*ssq2/ssqc)
     &		-log(u2*(1.0-u2**2)*u3*cu3)

	if(qdebug) then
		write(18,'(2i6,1pe16.7)') isweep,move,logratio
	end if

	logratio = max(-20.0,min(20.0,logratio))


	if(sdrand().lt.exp(-logratio)) then

c.. accept combine

		if(qdebug) then
		    write(18,'("combine accepted")')
		end if
	    naccc = naccc+1

c 	    next(j1) = next(j2)
c	    if(next(j1).ne.0) prev(next(j1)) = j1
	if(prev(j2).ne.0) then
		next(prev(j2)) = next(j2)
	else
		first = next(j2)
	end if
	if(next(j2).ne.0) prev(next(j2)) = prev(j2)
	    next(j2) = free
	    free = j2
	    prev(j2) = -1

	    leng(j1) = l1+l2
	    wt(j1) = wtc
	    mu(j1) = muc
	    ssq(j1) = ssqc
	    if(l1.ne.0) then
		inext(ilast) = st2
	    else
		start(j1) = st2
	    end if

	    k = k-1
	    if(l1*l2.eq.0) kemp = kemp-1
	    kpos = k-kemp

		if(.not.qkreg.and..not.qfix) write(8,'(2i3,i8)') k,kpos,isweep

	end if


	else if(move.eq.allocate) then

c	   write(0,*)'alloc'

c.. Gibbs move for allocation -------------------------------

	if(k.gt.1) then

	   call stdalloc(yv,n,wt,mu,ssq,ncmax,start,leng,next,pw,inext,
     &		first,delta,qprior)
	   kemp = 0
	   j = first
	   do while(j.ne.0)
	      if(leng(j).eq.0) kemp = kemp+1
	      j = next(j)
	   end do

	   if(qdebug) then
	      write(18,'("allocation move")')
	      j = first
	      do while(j.ne.0)
		 write(18,'(2i6)') j,leng(j)
		 j = next(j)
	      end do
	   end if
	   

	   
	end if
	
	else if(move.eq.weights) then

c	   write(0,*)'weights'

c.. Gibbs move for weights  ---------------------------------

	if(qdebug) then
		write(18,'("weights move")')
	end if
	wtsum = 0.0
	j = first
	do while(j.ne.0)
	    wt(j) = rgamma(delta+leng(j))
	    wtsum = wtsum+wt(j)
	    j = next(j)
	end do
	j = first
	do while(j.ne.0)
	    wt(j) = wt(j)/wtsum
	    if(qdebug) write(18,'(i4,f10.4)') j,wt(j)
	    j = next(j)
	end do


	else if(move.eq.parameters) then

c	   write(0,*)'params'

c.. Metropolis and/or Gibbs moves for component parameters --


	call gauss(mun,hwm)
	j = first
	do while(j.ne.0)
	if(qprior) then
	    if(qunif) then
		mun(j) = xi+unhw*(2.0*sdrand()-1.0)
	    else
		mun(j) = xi+mun(j)/sqrt(kappa)
	    end if
	else
	    ysum = 0.0
	    i = start(j)
	    do while(i.ne.0)
		ysum = ysum+yv(i)
		i = inext(i)
	    end do
	    if(qunif) then
		if(leng(j).eq.0) then
		    mun(j) = xi+unhw*(2.0*sdrand()-1.0)
		else
		    mun(j) = ysum/leng(j)+sqrt(ssq(j)/leng(j))*mun(j)
		end if
	    else
		con = 1.0/(leng(j)/ssq(j)+kappa)
		mun(j) = con*(ysum/ssq(j)+xi*kappa)+sqrt(con)*mun(j)
	    end if
	end if

	j = next(j)
	end do

c.. check order first
	jp = first
	jq = next(jp)
	do while(jq.ne.0)
	   if(mun(jp).ge.mun(jq)) go to 66
	   jp = jq
	   jq = next(jq)
	end do
	
c.. for unif option, have to check in range
	if(qunif) then
	   if(mun(first).lt.xi-unhw.or.mun(jp).gt.xi+unhw)
     &			go to 66
	end if

c.. if sp>1, calc acceptance prob, and reject if necessary
	if(sp.gt.1) then
	   sqrk = sqrt(kappa)
	   jp = first
	   temp = log(pnorm((mun(jp)-xi)*sqrk)/pnorm((mu(jp)-xi)*sqrk))
	   jq = next(jp)
	   do while(jq.ne.0)
	      temp = temp+log((pnorm((mun(jq)-xi)*sqrk)
     &		   -pnorm((mun(jp)-xi)*sqrk))
     &			/(pnorm((mu(jq)-xi)*sqrk)
     &			-pnorm((mu(jp)-xi)*sqrk)))
	      jp = jq
	      jq = next(jq)
	   end do
	   temp = temp+log((1.0-pnorm((mun(jp)-xi)*sqrk))
     &			/(1.0-pnorm((mu(jp)-xi)*sqrk)))
	   logratio = (sp-1)*temp
	   logratio = max(-20.0,min(0.0,logratio))
	   if(sdrand().gt.exp(logratio)) go to 66
	end if

c.. copy over if accepted
	j = first
	do while(j.ne.0)
	   mu(j) = mun(j)
	   j = next(j)
	end do
	

 66	if(qdebug) then
		write(18,'("means updated")')
		j = first
		do while(j.ne.0)
		    write(18,'(i4,f10.4)') j,mu(j)
		    j = next(j)
		end do
	end if

	j = first
	do while(j.ne.0)
	if(qprior) then
	    ssq(j) = beta/rgamma(alpha)
	else
	    ssd = 0.0
	    i = start(j)
	    do while(i.ne.0)
		ssd = ssd+(yv(i)-mu(j))**2
		i = inext(i)
	    end do
	    ssq(j) = (beta+0.5*ssd)/rgamma(alpha+0.5*leng(j))
	end if

	j = next(j)
	end do

	if(qdebug) then
		write(18,'("variances updated")')
		j = first
		do while(j.ne.0)
		    write(18,'(i4,f10.4)') j,ssq(j)
		    j = next(j)
		end do
	end if

	else if(move.eq.hyper) then

c	   write(0,*)'hyper'

c.. Gibbs move for hyperparameters --------------------------

	if(qbeta) then
		j = first
		sum = 0.0
		do while(j.ne.0)
			sum = sum+1.0/ssq(j)
			j = next(j)
		end do
		beta = rgamma(ggg+k*alpha)/(hhh+sum)
	end if
	if(qkappa) then
		j = first
		sum = 0.0
		do while(j.ne.0)
			sum = sum+(mu(j)-xi)**2
			j = next(j)
		end do
		kappa = rgamma(eee+0.5*k)/(fff+0.5*sum)
	end if


c.. Birth of an empty component

	else if(move.eq.birth) then

c	   write(0,*)'birth'

	    ntryb = ntryb+1

c.. compute logratio

	    klow = k
	    kemplow = kemp

	    wtstar = rbeta(1.0,real(klow))

	    ssqstar = beta/rgamma(alpha)

	    if(qunif) then
		mustar = xi+unhw*(2.0*sdrand()-1.0)
	    else
	    call gauss(mustar,1)
	    mustar = xi+mustar/sqrt(kappa)
	    end if

c.. p(k,w,z) terms

	    logratio = (lp(klow+1)-lp(klow))
     &		+(delta-1.0)*log(wtstar)-logbeta(delta,klow*delta)
     &		+(n+klow*(delta-1.0))*log(1.0-wtstar)

c.. p(theta) and proposal terms

	    logratio = logratio
     &		+log((klow+1)*d(klow+1)/(b(klow)*(kemplow+1)))
     &		-(klow-1)*log(1.0-wtstar)+logbeta(1.0,real(klow))

c.. Jacobian terms

	    logratio = logratio+(klow-1)*log(1.0-wtstar)

	    logratio = max(-20.0,min(20.0,logratio))

	    if(sdrand().lt.exp(logratio)) then

		scale = 1.0-wtstar
		j = first
		do while(j.ne.0)
		    wt(j) = scale*wt(j)
		    j = next(j)
		end do

		if(free.eq.0) stop 97
		jnew = free
		free = next(jnew)

		jprev = 0
		j = first
		do while(j.ne.0)
			if(mu(j).gt.mustar) go to 81
			jprev = j
			j = next(j)
		end do

81		next(jnew) = j
		prev(jnew) = jprev
		if(jprev.ne.0) then
			next(jprev) = jnew
		else
			first = jnew
		end if
		if(j.ne.0) prev(j) = jnew

		hwm = max(hwm,jnew)

		start(jnew) = 0
		leng(jnew) = 0
		mu(jnew) = mustar
		ssq(jnew) = ssqstar
		wt(jnew) = wtstar

		k = k+1
		kemp = kemp+1
		kpos = k-kemp

		if(.not.qkreg.and..not.qfix) write(8,'(2i3,i8)') k,kpos,isweep

		naccb = naccb+1

	    end if

c.. Death of an empty component

	else if(move.eq.death) then

c	   write(0,*)'death'

	    ntryd = ntryd+1

	    if(kemp.gt.0) then
80		jkill = 1+int(hwm*sdrand())
		if(prev(jkill).lt.0) go to 80
		if(leng(jkill).ne.0) go to 80

c.. compute logratio
		wtstar = wt(jkill)
		mustar = mu(jkill)
		ssqstar = ssq(jkill)
		klow = k-1
		kemplow = kemp-1

c.. p(k,w,z) terms

		logratio = (lp(klow+1)-lp(klow))
     &		+(delta-1.0)*log(wtstar)-logbeta(delta,klow*delta)
     &		+(n+klow*(delta-1.0))*log(1.0-wtstar)

c.. p(theta) and proposal terms

		logratio = logratio
     &		+log((klow+1)*d(klow+1)/(b(klow)*(kemplow+1)))
     &		-(klow-1)*log(1.0-wtstar)+logbeta(1.0,real(klow))

c.. Jacobian terms

		logratio = logratio+(klow-1)*log(1.0-wtstar)

		logratio = max(-20.0,min(20.0,logratio))

		if(sdrand().lt.exp(-logratio)) then

			if(prev(jkill).ne.0) then
				next(prev(jkill)) = next(jkill)
			else
				first = next(jkill)
			end if
			if(next(jkill).ne.0)
     &			    prev(next(jkill)) = prev(jkill)
			next(jkill) = free
			free = jkill
			prev(jkill) = -1

			scale = 1.0/(1.0-wtstar)
			j = first
			do while(j.ne.0)
			    wt(j) = scale*wt(j)
			    j = next(j)
			end do

			k = k-1
			kemp = kemp-1
			kpos = k-kemp

			if(.not.qkreg.and..not.qfix) write(8,'(2i3,i8)') k,kpos,isweep

			naccd = naccd+1

		end if

	   end if


c.. Trap for illegal move

	else

		stop 888

	end if



c.. end of 'istep' loop:
	end do

c.. logging stage

	kpos = k-kemp

	if(qkreg.and..not.qfix) then
	    if(mod(isweep,nskdel).eq.0) then
		write(8,'(2i3,i8)') k,kpos,isweep
	    end if
	end if

	if(isweep.gt.0) then

	count(k) = count(k)+1
	countpos(kpos) = countpos(kpos)+1

c.. updates complete: record information

	if(qrkpos) then
		krep = kpos
	else
		krep = k
	end if
	if(krep.le.ncmax2) then
		j = first
		ij = 1
		do while(j.ne.0)
		    if((.not.qrkpos).or.leng(j).gt.0) then
			wtav(krep,ij) = wtav(krep,ij)+wt(j)
			muav(krep,ij) = muav(krep,ij)+mu(j)
			sigav(krep,ij) = sigav(krep,ij)+sqrt(ssq(j))
			avn(krep,ij) = avn(krep,ij)+leng(j)
			ij = ij+1
		    end if
		    j = next(j)
		end do
	end if

	if((.not.qprior).and.qfull) then

	    dev = 0.0
	    devc = 0.0
	    do i = 1,n
	       dvdy = 1.0
		dens = 0.0
		dens2 = 0.0
		j = first
		do while(j.ne.0)
		    temp = rr2pi*
     &			exp(max(-20.0,-0.5*(yv(i)-mu(j))**2/ssq(j)))/
     &			sqrt(ssq(j))
		    dens = dens+wt(j)*temp
		    dens2 = dens2+((leng(j)+delta)/(n+k*delta))*temp
		    j = next(j)
		end do
		dev = dev+log(dens*dvdy)
		devc = devc+log(dens2*dvdy)
	    end do
	    dev = -2.0*dev
	    devc = -2.0*devc
c	    devc = dev+devcorrection

	    avdev(krep) = avdev(krep)+dev
	    avdevc(krep) = avdevc(krep)+devc

	end if

	if(mod(isweep,nspace).eq.0) then

	ent = entropy(n,leng,ncmax,first,next)

	if(qbeta.or.qkappa) write(10,*) beta,kappa,k,kpos

	if(qfull) write(9,'(2i4,2f8.2)') k,kpos,dev,devc

	if(krep.le.ncmax2.and.qfull.and.(qpwms.or.qpalloc)) then
	    if(.not.fileopen(krep)) then
		istd2 = 6-int(log10(0.5+krep))
		write(num2,'(i6)') krep
		if(qpwms) then
		open(25+krep,file=
     &			prfx(1:npf)//".wms."//num2(istd2:6),
     &                  status='unknown')
		end if
		if(qpalloc) then
		open(25+ncmax+krep,file=
     &			prfx(1:npf)//".z."//num2(istd2:6),
     &                  status='unknown')
		end if
		fileopen(krep) = .true.
	    end if
	    if(qpwms) then
	        j = first
		do while(j.ne.0)
		    if((.not.qrkpos).or.leng(j).gt.0) then
			write(25+krep,'(3f12.5)') wt(j),mu(j),sqrt(ssq(j))
		    end if
		    j = next(j)
		end do
	    end if
	    if(qpalloc) then
		j = first
		ij = 1
		do while(j.ne.0)
			i = start(j)
			do while(i.ne.0)
				z(i) = ij
				i = inext(i)
			end do
			j = next(j)
			ij = ij+1
		end do
		write(25+ncmax+krep,'(20i3)') (z(i),i=1,n)
	    end if
	end if

	end if

	nkemp = nkemp+1
	avkemp = avkemp+kemp
	if(kemp.gt.0) ppkemp = ppkemp+1.0
	if(kemp.le.9) nfkemp(kemp) = nfkemp(kemp)+1

	if(mod(isweep,nspace).eq.0) then
		write(17,*) k,kpos,ent
	end if

	if(mod(isweep,max(1,nsweep/nsamp)).eq.0) then

	write(6,'(2i3)') k,kpos
	j = first
	do while(j.ne.0)
	    stdev = sqrt(ssq(j))
	    write(6,'(3f12.5,i4)') wt(j),mu(j),stdev,leng(j)
	    j = next(j)
	end do
	end if

c.. accumulate mean densities

	if(qfull) then
	    krec = min(k,ncd)
	    countden(krec) = countden(krec)+1
	    do iyd = 1,ngrid
	       yvg = yd0+iyd*ydinc
	       dvdy = 1.0
		pwsum = 0.0
		j = first
		ij = 1
		do while(j.ne.0)
		    pw(ij) = rr2pi*wt(j)*
     &			exp(max(-20.0,-0.5*(yvg-mu(j))**2/
     &			ssq(j)))/sqrt(ssq(j))
		    den(krec,iyd) = den(krec,iyd)+pw(ij)*dvdy
		    pwsum = pwsum+pw(ij)
		    j = next(j)
		    ij = ij+1
		end do
		if(qpclass) then
		    if(k.ge.k1.and.k.le.k2) then
		    ko = ((k+k1-3)*(k-k1))/2
		    do ij = 1,k-1
			pclass(ko+ij,iyd) = pclass(ko+ij,iyd)
     &				+pw(ij)/pwsum
		    end do
		    end if
		end if
	    end do

		if(qpclass) then
		if(k.ge.k1.and.k.le.k2) then
		    ko = ((k+k1-3)*(k-k1))/2
		    j = first
		    ij = 1
		    do while(ij.lt.k)
			i = start(j)
			do while(i.ne.0)
			    pz(ko+ij,i) = pz(ko+ij,i)+1
			    i = inext(i)
			end do
			j = next(j)
			ij = ij+1
		    end do
		end if
		end if
	end if

	if(qstop) stop 99

	end if

c.. end of main loop

	if(mod(isweep,100).eq.0) then
		hwm = 0
		j = first
		do while(j.ne.0)
			hwm = max(hwm,j)
			j = next(j)
		end do
	end if

	if(qdebug.and.mod(isweep,nsweep/50).eq.0) write(0,*) isweep

	end do

	write(0,*)

	write(7,'("splits:  ",i6," out of",i7)') naccs,ntrys
	write(7,'("split rej r:  ",i6)') nrejr
	write(7,'("combines:",i6," out of",i7)') naccc,ntryc
	if(qempty) then
		write(7,'("births:  ",i6," out of",i7)') naccb,ntryb
		write(7,'("deaths:  ",i6," out of",i7)') naccd,ntryd
	end if

c.. output posterior probs for k, and Bayes factors

	kbase = 0
	do k = 1,ncmax
		pw(k) = real(count(k))/nsweep
		if(kbase.eq.0.and.count(k).ne.0) kbase = k
	end do
	write(12,'(5f16.6)') pw

	do k = 1,ncmax
	    bf(k) = (pw(k)/pw(kbase))/exp(lp(k)-lp(kbase))
	end do

	write(12,'(1x)')
	write(12,'(5f16.6)') bf

c.. write prior probs p(k)

	do k = 1,ncmax
	    pw(k) = exp(lp(k))
	end do
	write(12,'(1x)')
	write(12,'(5f16.6)') pw

c.. write posterior probs for number of nonempty components

	do k = 1,ncmax
		   pw(k) = real(countpos(k))/nsweep
	end do
	write(12,'(1x)')
	write(12,'(5f16.6)') pw

	avkemp = avkemp/nkemp
	ppkemp = ppkemp/nkemp
	write(12,'(2f10.4)') avkemp,ppkemp
	write(12,'(10i8)') nfkemp

	if(qfull) then

	    avdevall = 0.0
	    avdevcall = 0.0
	    do k = 1,ncmax
		avdevall = avdevall+avdev(k)
		avdevcall = avdevcall+avdevc(k)
	    end do
	    avdevall = avdevall/nsweep
	    avdevcall = avdevcall/nsweep

	    if(qrkpos) then
		do k = 1,ncmax
		    avdev(k) = avdev(k)/max(1,countpos(k))
		    avdevc(k) = avdevc(k)/max(1,countpos(k))
		end do
	    else
		do k = 1,ncmax
		    avdev(k) = avdev(k)/max(1,count(k))
		    avdevc(k) = avdevc(k)/max(1,count(k))
		end do
	    end if
	    write(12,'(1x)')
	    write(12,'(5f16.6)') avdev
	    write(12,'(1x)')
	    write(12,'(5f16.6)') avdevc
	end if

c.. output posterior expectations of parameters

	do k = 1,ncmax2
	if(qrkpos) then
		ckrep = countpos(k)
	else
		ckrep = count(k)
	end if
	write(11,'(i3,i8)') k,ckrep
	ckrep = max(1,ckrep)
	do ij = 1,k
		wtav(k,ij) = wtav(k,ij)/ckrep
		muav(k,ij) = muav(k,ij)/ckrep
		sigav(k,ij) = sigav(k,ij)/ckrep
		write(11,'(3f12.5)') wtav(k,ij),muav(k,ij),sigav(k,ij)
	end do
	end do

c.. output average group sizes

	if(qfull) then
	    do k = 1,ncmax2
		if(qrkpos) then
		    ckrep = countpos(k)
		else
		    ckrep = count(k)
		end if
		write(16,'(2i8)') k,ckrep
		ckrep = max(1,ckrep)
		do ij = 1,k
		    avn(k,ij) = avn(k,ij)/ckrep
		end do
		write(16,'(10f8.2)') (avn(k,ij),ij=1,k)
	    end do
	end if
c.. output densities

	if(qfull) then
	    do k = 1,ncd-1
	    countden(ncd) = countden(ncd)+countden(k)
	    do iyd = 1,ngrid
	    den(ncd,iyd) = den(ncd,iyd)+den(k,iyd)
	    end do
	    end do
	    if(qfull)then
		write(13,*) (yd0+iyd*ydinc,iyd=1,ngrid)
	    end if
	    do k = 1,ncd
	    do iyd = 1,ngrid
	    den(k,iyd) = den(k,iyd)/max(countden(k),1)
	    end do
	    if(qfull)then
		write(13,*) (den(k,iyd),iyd=1,ngrid)
	    end if
	    end do

	    write(12,'("     D-hat:",g13.5)') avdevall
	    write(12,'(" new D-hat:",g13.5)') avdevcall
	    do k = 1,ncd
		ggh(k) = 0.0
		if(countden(k).gt.0) then
		  do i = 1,n
		    iyd = max(1,min(ngrid,int(0.5+(y(i)-yd0)/ydinc)))
		    ggh(k) = ggh(k)+log(den(k,iyd))
		  end do
		  ggh(k) = -2.0*ggh(k)
		end if
	    end do
	    write(12,'("  G(g-hat):",g13.5)') ggh(ncd)
	    write(12,'("G(g-hat_k):",3g13.5/(11x,3g13.5))')
     &		(ggh(k),k=1,ncd-1)

	end if

c.. output predictive classification

	if(qfull.and.qpclass) then

	write(14,*) (yd0+iyd*ydinc,iyd=1,ngrid)
	do k = k1,k2
	    ko = ((k+k1-3)*(k-k1))/2
	    write(14,'(i4)') k
	    do ij = 1,k-1
		do iyd = 1,ngrid
		    pclass(ko+ij,iyd) = pclass(ko+ij,iyd)/max(1,count(k))
		end do
		if(ij.ne.1) then
		    do iyd = 1,ngrid
			pclass(ko+ij,iyd) = pclass(ko+ij,iyd)+
     &				pclass(ko+ij-1,iyd)
		    end do
		end if
	    write(14,'(10f8.4)') (pclass(ko+ij,iyd),iyd=1,ngrid)
	    end do
	end do

c.. output within-sample classification

	do k = k1,k2
	    ko = ((k+k1-3)*(k-k1))/2
	    write(15,'(i4)') k
	    do ij = 1,k-1
		do i = 1,n
		    pz(ko+ij,i) = pz(ko+ij,i)/max(1,count(k))
		end do
		if(ij.ne.1) then
		    do i = 1,n
			pz(ko+ij,i) = pz(ko+ij,i)+pz(ko+ij-1,i)
		    end do
		end if
		write(15,'(10f8.4)') (pz(ko+ij,i),i=1,n)
	    end do
	end do

	end if

	stop
	end

c---------------------------------------------------

	subroutine stdalloc(y,n,wt,mu,ssq,ncmax,start,leng,next,pw,inext,
     &		first,delta,qprior)

	integer first,start
	real mu
	logical qprior
	dimension start(ncmax),leng(ncmax),inext(n),next(ncmax),
     &		pw(ncmax),y(n),wt(ncmax),mu(ncmax),ssq(ncmax)

	j = first
	do while(j.ne.0)
 	    start(j) = 0
	    leng(j) = 0
	    j = next(j)
	end do

	do i = 1,n

	pwsum = 0.0

	j = first
	do while(j.ne.0)
	if(qprior) then
	    pw(j) = wt(j)
	else
	    pw(j) = wt(j)
     &		*exp(max(-20.0,-0.5*(y(i)-mu(j))**2/ssq(j)))/sqrt(ssq(j))
	end if
	    pwsum = pwsum+pw(j)
	    j = next(j)
	end do
	usd = sdrand()
	u = pwsum*usd

	j = first
	do while(j.ne.0)
	    u = u-pw(j)
	    if(u.lt.0.0) go to 43
	    j = next(j)
	end do

	j = first

43	inext(i) = start(j)
	start(j) = i
	leng(j) = leng(j)+1

	end do

	return

	end

c---------------------------------------------------

	subroutine reorder(mu,ncmax,next,prev,first)

c.. check ordering of mu's and correct if necessary

	integer prev,first
	real mu

	dimension mu(ncmax),prev(ncmax),next(ncmax)

	jp = first
73	jnext = next(jp)
	jq = prev(jp)
74	if(jq.eq.0) go to 75
	if(mu(jq).lt.mu(jp)) go to 75
	jq = prev(jq)
	go to 74
75	if(prev(jp).ne.jq) then
		if(prev(jp).le.0) stop 90
c.. remove jp from list
		next(prev(jp)) = jnext
		if(jnext.ne.0) prev(jnext) = prev(jp)
c.. re-insert jp, after jq
		prev(jp) = jq
		if(jq.ne.0) then
			next(jp) = next(jq)
			if(next(jq).le.0) stop 92
			prev(next(jq)) = jp
			next(jq) = jp
		else
			next(jp) = first
			prev(first) = jp
			first = jp
		end if
	end if
	jp = jnext
	if(jp.ne.0) go to 73

	return

	end

c---------------------------------------------------

      subroutine sort(na,llo,lhi,q,nq)
      dimension na(lhi),q(nq)
c
c    sorts array na(l), l = llo to lhi
c    according to values of q(na(l))
c    so if llo = 1 and na(l) = l on entry,
c    obtain na(l) = index of l-th smallest q on exit.
c
      mesh = lhi-llo+1
    1 mesh = (mesh+1)/3
      lst = llo+mesh
      do 3 l1 = lst,lhi
      nal = na(l1)
      ql = q(nal)
      lnow = l1
      do 2 l2 = lst,l1,mesh
      lnext = lnow-mesh
      nan = na(lnext)
      if(ql.ge.q(nan)) go to 3
      na(lnow) = nan
    2 lnow = lnext
    3 na(lnow) = nal
      if(mesh.gt.1) go to 1
      return
      end

c---------------------------------------------------

	function entropy(n,leng,ncmax,first,next)

c.. calculates ent = -\sum_j (n_j/n) \log (n_j/n)

	integer first,leng(ncmax),next(ncmax)

	entropy = log(real(n))

	j = first
	do while(j.ne.0)
		if(leng(j).ne.0) then
		    temp = real(leng(j))
		    entropy = entropy-(temp/n)*log(temp)
		end if
		j = next(j)
	end do

	return

	end

c---------------------------------------------------

	subroutine rmu(mu,y,xi,kappa,alpha,beta)

c		mu = rv from f propto
c                   exp(-kappa(mu-xi)^2/2)/(beta+.5*(yi-mu)^2)^(alpha+.5)

	real mu,kappa

1	xx = sdrand()
	yy = xx*y+0.2*(sdrand()-0.5)
	mu = yy/xx
	if(xx**2.gt.exp(-0.5*kappa*(mu-xi)**2)/
     &	    (1.0+(0.5/beta)*(y-mu)**2)**(alpha+0.5)) go to 1

	return

	end

c---------------------------------------------------

	function rbeta(a1,a2)
	x = rgamma(a1)
	y = rgamma(a2)
	rbeta = x/(x+y)
	return
	end

c---------------------------------------------------

	real function logbeta(x1,x2)
	real*8 dlgama
	logbeta = sngl(dlgama(dble(x1))+dlgama(dble(x2))
     &		-dlgama(dble(x1+x2)))
	return
	end
