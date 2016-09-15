	subroutine gauss(z,n)

c..    generates an n-vector of i.i.d. N(0,1) r.v.'s z
c      by Box-Mueller method

	real z(n)

	n1 = n-1
	do 1 i = 1,n1,2
	u = sqrt(-2.0*alog(sdrand(0)))
	v = 6.2831853*sdrand(0)
	z(i) = u*sin(v)
1	z(i+1) = u*cos(v)
	if(mod(n,2).eq.0) return
	z(n) = sqrt(-2.0*alog(sdrand(0)))*sin(6.2831853*sdrand(0))

	return
	end
