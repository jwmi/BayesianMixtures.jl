	function rgamma(s)
	data e/2.71828182/

c.. Returns one value from the gamma distribution
c   GAMMA(s,1), i.e. with pdf = y^(s-1) exp(-y) / GAMMA(s)

	if(s.lt.1) then
		b = (s+e)/e
		c1 = 1.0/s
1		bu = b*sdrand()
		if(bu.le.1.0) then
c.. this mod is to prevent underflow when s is small
		    rgamma = exp(max(-30.0,c1*log(bu)))
		    if(sdrand().lt.exp(-rgamma)) return
		    go to 1
		else
		    rgamma = -log((b-bu)/s)
		    if(sdrand().lt.rgamma**(s-1)) return
		    go to 1
		end if
	else if(s.eq.1) then
		rgamma = -log(sdrand())
		return
	else
		c1 = s-1.0
		c2 = (s-1.0/(6.0*s))/c1
		c3 = 2.0/c1
		c4 = c3+2.0
		c5 = 1.0/sqrt(s)
2		u1 = sdrand()
		u2 = sdrand()
		if(s.gt.2.5) u1 = u2+c5*(1.0-1.86*u1)
		if(u1.le.0.0.or.u1.ge.1.0) go to 2
		w = c2*u2/u1
		if(c3*u1+w+1.0/w.le.c4) go to 3
		if(c3*log(u1)-log(w)+w.ge.1.0) go to 2
3		rgamma = c1*w
		return
	end if

	end
