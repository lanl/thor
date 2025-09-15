#!/usr/bin/env python3
# %%

import numpy
import math

class BaseQuadRule:
	## attributes
	id = ""
	nq = 0
	wq = []
	sq = []
	## methods
	def get_id(self):
		return self.id
	def get_nq(self):
		return self.nq
	def get_wq(self,iq):
		if ( iq<self.nq and self.nq>0 ): 
			return self.wq[iq]
	def get_sq(self,iq):
		if ( iq<self.nq and self.nq>0 ): 
			return self.sq[iq]
	def get_quad(self,iq,a,b):
		if ( iq<self.nq and self.nq>0 ):
			wg = (b-a)*self.wq[iq]/2.
			xg = a + (b-a)*(1.+self.sq[iq])/2.
			return [wg,xg]
	def get_wqvec(self):
		return numpy.array(self.wq)
	def get_weight_and_grid(self, a, b, nc):
		dx = (b - a)/nc
		nq = self.get_nq()
		W = []
		X = []
		for i0 in range(nc):  # loop on cells
			x0 = a + i0 * dx
			x1 = x0 + dx
			for iq in range(nq):  # accumulate on quad nodes
				[wq, xq] = self.get_quad(iq, x0, x1)
				W.append(wq)
				X.append(xq)
		return numpy.array(W), numpy.array(X) 
				
## mid point rule
class MidPointRule(BaseQuadRule):
		id = "MidPointRule"
		nq = 1
		wq = [ 2. ] ## multiplied by 2
		sq = [ 0. ]

## trapezoidal rule
class Trapezoidal(BaseQuadRule):
		id = "Trapezoidal"
		nq = 2
		wq = [  1.0, 1.0 ] ## multiplied by 2
		sq = [ -1.0, 1.0 ]

## Simpson rule
class Simpson(BaseQuadRule):
		id = "Simpson"
		nq = 3
		wq = [  2./6., 8./6., 2./6. ] ## multiplied by 2
		sq = [ -1.0,   0.0,   1.0   ]
		
## Gauss2
class Gauss2(BaseQuadRule):
		id = "Gauss2"
		nq = 2
		wq = [ 1., 1. ]
		sq = [ -1./math.sqrt(3.), +1./math.sqrt(3) ]

## Gauss3
class Gauss3(BaseQuadRule):
		id = "Gauss3"
		nq = 3
		wq = [ 5./9., 8./9., 5./9. ]
		sq = [ -0.774596669241483, +0.000000000000000, +0.774596669241483 ]

## Gauss4
class Gauss4(BaseQuadRule):
		id = "Gauss4"
		nq = 4
		wq = [  0.347854845137454,  0.652145154862546,  0.652145154862546,  0.347854845137454 ]
		sq = [ -0.861136311594053, -0.339981043584856, +0.339981043584856, +0.861136311594053 ]

## Gauss8
class Gauss8(BaseQuadRule):
		id = "Gauss8"
		nq = 8
		wq = [  0.1012285362903760,  0.2223810344533743,  0.3137066458778873,  0.3626837833783622,\
						0.3626837833783619,  0.3137066458778869,  0.2223810344533742,  0.1012285362903759 ]
		sq = [ -0.9602898564975365, -0.7966664774136270, -0.5255324099163290, -0.1834346424956498,\
					 +0.1834346424956496, +0.5255324099163292, +0.7966664774136268, +0.9602898564975364 ]

## Gauss16
class Gauss16(BaseQuadRule):
		id = "Gauss16"
		nq = 16
		wq = [  0.02715245941175414,  0.06225352393864799,  0.09515851168249308,  0.12462897125553410,\
						0.14959598881657670,  0.16915651939500240,  0.18260341504492330,  0.18945061045506820,\
						0.18945061045506820,  0.18260341504492330,  0.16915651939500240,  0.14959598881657680,\
						0.12462897125553410,  0.09515851168249308,  0.06225352393864799,  0.02715245941175415 ]
		sq = [ -0.98940093499164990, -0.94457502307323260, -0.86563120238783190, -0.75540440835500310,\
					 -0.61787624440264380, -0.45801677765722740, -0.28160355077925890, -0.09501250983763744,\
					 +0.09501250983763744, +0.28160355077925890, +0.45801677765722740, +0.61787624440264380,\
					 +0.75540440835500310, +0.86563120238783190, +0.94457502307323260, +0.98940093499164990 ]
		
## FUNCTIONS

## check the quadrule
def prt_quadrule_1D( qd ):
		print(qd.get_id())
		nq = qd.get_nq()
		for iq in range(nq):
				wq = qd.get_wq(iq)
				sq = qd.get_sq(iq)
				print("iq = %d  wq = %24.16e  sq = %24.16e" %(iq,wq,sq))
		# ---
		a  = 0
		b  = 1
		np = 1
		sum_q = compute_intg_poly(np,qd,a,b)
		print("sum_q = %24.16e" %(sum_q))
		print("---------------------------------------------------------------------------------------")

## compute 1D integral on interval [a,b]
def compute_intg_poly(np,qd,a,b):
		ng = qd.get_nq()
		sum_q = 0.
		fq = lambda xq : (np+1)*xq**np
		for iq in range(ng):
				[wq,xq] = qd.get_quad(iq,a,b)
				sum_q += wq * fq(xq)
		return sum_q
## -->> end of compute_intg_poly

## evaluate approximation error and convergence rate
def evaluate_intg_convg_rate(qd,nd,a,b,np_min=0,np_max=0):
		ndx = [10, 20, 30, 40] #!number of intervals ???
		val = [ 0., 0., 0., 0., ]
		err = [ 0., 0., 0., 0., ]
		if ( np_max==0 ): np_max = 10
		for ip in range(np_min,np_max+1,1):
				i = 0
				print()
				print("polynomial degree = ",ip)
				for ni in ndx:
						nc = []
						for j in range(nd):
								nc.append(ni)
						## ---------------
						val[i] = compute_intg_poly_1D(ip,qd,nd,a,b,nc)
						err[i] = abs( 1. - val[i] ) + 1.e-20
						if ( ni>10 ):
								rate = math.log( err[i]/err[i-1] ) / math.log( ndx[i-1]/ndx[i] )
								print("ni = %d val = %24.16e -->> err = %14.7e  rate = %14.7e" %(ni,val[i],err[i],rate))
						else:
								print("ni = %d val = %24.16e -->> err = %14.7e  rate = --- " %(ni,val[i],err[i]))
						## --> end of if ( ni>0 )...
						i += 1
				## -->> end of for ni in ndx:
		## -->> end of for ip in range(np):
## -->> end of evaluate_intg_convg_rate_1D


## compute 1D integral on mesh of [a,b]
def compute_intg_poly_1D(np,qd,nd,a,b,nc):
		## function to be integrated
		fq = lambda xq : (np+1)*xq**np
		## ----
		dx = ( b[0] - a[0]) / float(nc[0])
		## ----
		nq = qd.get_nq()
		sum_q = 0.
		W = []
		X = []
		for i0 in range(nc[0]):       ## loop on cells
				x0 = i0 * dx
				x1 = x0 + dx
				for iq in range(nq):   ## accumulate on quad nodes
					[wq,xq] = qd.get_quad(iq,x0,x1)
					W.append(wq)
					X.append(xq)
					sum_q += wq * fq(xq)
				# end of for iq...
		#end of -- for ic...
		return sum_q
## -->> end of compute_intg_poly_1D

## evaluate integrals in multi-dimensional setting
def compute_intg_poly_2D(np,qd,nd,a,b,nc):
		# print("np = ",np)
		# print("nd = ",nd)
		# print("a  = ",a)
		# print("b  = ",b)
		# print("nc = ",nc)
		## function to be integrated
		func = lambda x0, x1 : float(np+1)*( x0**np + x1**np )/float(nd)
		## ----
		dx = [ 0., 0. ]
		xb = [ 0., 0. ]
		xe = [ 0., 0. ]
		## ----
		for i in range(nd):
				dx[i] = ( b[i] - a[i] ) / float(nc[i])
		## ----
		wq = []
		fq = []
		nq = qd.get_nq()
		iq = 0
		for i0 in range(nc[0]):       ## loop on 1D partition along direction i1
				xb[0] = i0 * dx[0]
				xe[0] = xb[0] + dx[0]
				for i1 in range(nc[1]):       ## loop on 1D partition along direction i2
						xb[1] = i1 * dx[1]
						xe[1] = xb[1] + dx[1]
						for iq0 in range(nq):       ## loop on quad nodes along direction i0
								for iq1 in range(nq):   ## loop on quad nodes along direction i1
										[wq0,xq0] = qd.get_quad(iq0,xb[0],xe[0])
										[wq1,xq1] = qd.get_quad(iq1,xb[1],xe[1])
										wq.append( wq0 * wq1 )
										fq.append( func( xq0, xq1 ) )
										iq += 1
								# end of -->> for iq1...
						# end of -->> for iq0...
				#end of -->> for i1...
		#end of -->> for i0...

		## set total number of quad terms in 2D
		nq_tot = iq
		
		## accumulate quad rule evaluations
		sum_q = 0.
		for iq in range(nq_tot):
				sum_q += wq[iq] * fq[iq]
		## end of for iq...
		
		return sum_q
## -->> end of compute_intg_poly_2D

## evaluate integrals in multi-dimensional setting
def compute_intg_poly_3D(np,qd,nd,a,b,nc):
		## function to be integrated
		func = lambda x0, x1, x2 : float(np+1)*( x0**np + x1**np + x2**np)/float(nd)
		## ----
		dx = [ 0., 0., 0. ]
		xb = [ 0., 0., 0. ]
		xe = [ 0., 0., 0. ]
		## ----
		for i in range(nd):
				dx[i] = ( b[i] - a[i] ) / float(nc[i])
		## ----
		wq = []
		fq = []
		nq = qd.get_nq()
		iq = 0
		for i0 in range(nc[0]):       ## loop on 1D partition along direction i1
				xb[0] = i0 * dx[0]
				xe[0] = xb[0] + dx[0]
				for i1 in range(nc[1]):       ## loop on 1D partition along direction i2
						xb[1] = i1 * dx[1]
						xe[1] = xb[1] + dx[1]
						for i2 in range(nc[2]):       ## loop on 1D partition along direction i2
								xb[2] = i2 * dx[2]
								xe[2] = xb[2] + dx[2]
								for iq0 in range(nq):       ## loop on quad nodes along direction i0
										for iq1 in range(nq):   ## loop on quad nodes along direction i1
												for iq2 in range(nq):   ## loop on quad nodes along direction i2
														[wq0,xq0] = qd.get_quad(iq0,xb[0],xe[0])
														[wq1,xq1] = qd.get_quad(iq1,xb[1],xe[1])
														[wq2,xq2] = qd.get_quad(iq2,xb[2],xe[2])
														wq.append( wq0 * wq1 * wq2 )
														fq.append( func( xq0, xq1, xq2 ) )
														iq += 1
												# end of -->> for iq2...
								# end of -->> for iq1...
						# end of -->> for iq0...
				#end of -->> for i1...
		#end of -->> for i0...

		## set total number of quad terms in 2D
		nq_tot = iq
		
		## accumulate quad rule evaluations
		sum_q = 0.
		for iq in range(nq_tot):
				sum_q += wq[iq] * fq[iq]
		## end of for iq...
		
		return sum_q
## -->> end of compute_intg_poly_3D


	# %% 
if __name__ == '__main__':
	# test compute_intg_poly_1D
	# qd = Trapezoidal() #quadrature rule
	# nd = 1 # number of dimension
	# a = [0]
	# b = [1]
	# nc = [10] # number of intervals
	# np = 0 # polynomial degree
	# s, W, X = compute_intg_poly_1D(np,qd,nd,a,b,nc)
	# print(s)
	# print(W)
	# print(X)
	
	
			
	# %% 
	# # just an abbreviation
	T = True

	## ------------------------------------------
	## MAIN - test 1D numerical integration
	## ------------------------------------------

	##              Trap   Simp   MidPt  GL-2   GL-3   GL-4   GL-8   GL-16
	##              0      1      2      3      4      5      6      7
	run_flag_1D = [ False, False, False, False, False,  False,  False, T]

	nd = 1
	a  = [ 0. ]
	b  = [ 1. ]

	## ------------------------------------------
	## NEWTON-COTES formulas (including extremes)
	## ------------------------------------------
	# %%
	## trapezoidal rule exact on polynomials of degree 1
	if ( run_flag_1D[0] ):
			qd = Trapezoidal()
			prt_quadrule_1D(qd)
			evaluate_intg_convg_rate(qd,nd,a,b)

	## Simpson rule exact on polynomials of degree 3
	if ( run_flag_1D[1] ):
			qd = Simpson()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b)

	## ------------------------------------------------
	## GAUSS-LEGENDRE formulas (not including extremes)
	## ------------------------------------------------

	## mid-point rule exact on polynomials of degree 1
	if ( run_flag_1D[2] ):
			qd = MidPointRule()
			prt_quadrule_1D(qd)
			evaluate_intg_convg_rate(qd,nd,a,b)

	## 2-pts Gauss-Legendre rule exact on polynomials of degree 3
	if ( run_flag_1D[3] ):
			qd = Gauss2()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b)

	## 3-pts Gauss-Legendre rule exact on polynomials of degree 5
	if ( run_flag_1D[4] ):
			qd = Gauss3()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,2,7)

	## 4-pts Gauss-Legendre rule exact on polynomials of degree 7
	if ( run_flag_1D[5] ):
			qd = Gauss4()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,4,9)

	## 8-pts Gauss-Legendre rule exact on polynomials of degree 15
	if ( run_flag_1D[6] ):
			qd = Gauss8()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,12,17)

	## 16-pts Gauss-Legendre rule exact on polynomials of degree 31
	if ( run_flag_1D[7] ):
			qd = Gauss16()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,28,33)


	## ------------------------------------------
	## MAIN - test 2D numerical integration
	## ------------------------------------------

	##              Trap   Simp   MidPt  GL-2   GL-3   GL-4   GL-8   GL-16
	##              0      1      2      3      4      5      6      7
	run_flag_2D = [ False, False, False, False, False,  False,  False, False]

	nd = 2
	a  = [ 0., 0. ]
	b  = [ 1., 1. ]

	## ------------------------------------------
	## NEWTON-COTES formulas (including extremes)
	## ------------------------------------------

	## trapezoidal rule exact on polynomials of degree 1
	if ( run_flag_2D[0] ):
			qd = Trapezoidal()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate_2D(qd)

	## Simpson rule exact on polynomials of degree 3
	if ( run_flag_2D[1] ):
			qd = Simpson()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate_2D(qd)

	## ------------------------------------------------
	## GAUSS-LEGENDRE formulas (not including extremes)
	## ------------------------------------------------

	## mid-point rule exact on polynomials of degree 1
	if ( run_flag_2D[2] ):
			qd = MidPointRule()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate_2D(qd)
			evaluate_intg_convg_rate(qd,nd,a,b)
			
	## 2-pts Gauss-Legendre rule exact on polynomials of degree 3
	if ( run_flag_2D[3] ):
			qd = Gauss2()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,0,5)
			
	## 3-pts Gauss-Legendre rule exact on polynomials of degree 5
	if ( run_flag_2D[4] ):
			qd = Gauss3()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,2,7)

	## 4-pts Gauss-Legendre rule exact on polynomials of degree 7
	if ( run_flag_2D[5] ):
			qd = Gauss4()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,5,11)

	## 8-pts Gauss-Legendre rule exact on polynomials of degree 15
	if ( run_flag_2D[6] ):
			qd = Gauss8()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,14,19)

	## 16-pts Gauss-Legendre rule exact on polynomials of degree 31
	if ( run_flag_2D[7] ):
			qd = Gauss16()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,28,33)


	## ------------------------------------------
	## MAIN - test 3D numerical integration
	## ------------------------------------------

	##              Trap   Simp   MidPt  GL-2   GL-3   GL-4   GL-8   GL-16
	##              0      1      2      3      4      5      6      7
	run_flag_3D = [ False, False, False, False, False,  False,  False, False]

	nd = 3
	a  = [ 0., 0., 0. ]
	b  = [ 1., 1., 1. ]

	## ------------------------------------------
	## NEWTON-COTES formulas (including extremes)
	## ------------------------------------------

	## trapezoidal rule exact on polynomials of degree 1
	if ( run_flag_3D[0] ):
			qd = Trapezoidal()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b)

	## Simpson rule exact on polynomials of degree 3
	if ( run_flag_3D[1] ):
			qd = Simpson()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,np_min=0,np_max=0)

	## ------------------------------------------------
	## GAUSS-LEGENDRE formulas (not including extremes)
	## ------------------------------------------------

	## mid-point rule exact on polynomials of degree 1
	if ( run_flag_3D[2] ):
			qd = MidPointRule()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b)

	## 2-pts Gauss-Legendre rule exact on polynomials of degree 3
	if ( run_flag_3D[3] ):
			qd = Gauss2()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,0,5)
			
	## 3-pts Gauss-Legendre rule exact on polynomials of degree 5
	if ( run_flag_3D[4] ):
			qd = Gauss3()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,2,7)
			
	## 4-pts Gauss-Legendre rule exact on polynomials of degree 7
	if ( run_flag_3D[5] ):
			qd = Gauss4()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,5,10)
			
	## 8-pts Gauss-Legendre rule exact on polynomials of degree 15
	if ( run_flag_3D[6] ):
			qd = Gauss8()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,12,17)
			
	## 16-pts Gauss-Legendre rule exact on polynomials of degree 31
	if ( run_flag_3D[7] ):
			qd = Gauss16()
			prt_quadrule_1D( qd )
			evaluate_intg_convg_rate(qd,nd,a,b,28,33)
			
