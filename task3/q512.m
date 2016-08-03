%% Question 3.5
% The code for simulating Examples 3.4 and 3.5 is listed below. Study the 
% code and try to understand the details. 
%%

% 1.	 DISCRETE SYSTEM reg
% 2.	 "Indirect Self-Tuning Regulator based on the model
% 3.	 " H(q)=(b0*q+b1)/(q^2+a1*q+a2)
% 4.	 "using standard RLS estimation and pole placement design
% 5.     "Polynomial B is canceled if cancel>0.5
% 6.	 INPUT ysp y "set point and process output
% 7.	 OUTPUT u "control variable
% 8.	 STATE ysp1 y1 u1 v1 "controller states
% 9.	 STATE th1 th2 th3 th4 "parameter estimates
% 10.	 STATE f1 f2 f3 f4 "regression variables
% 11.	 STATE p11 p12 p13 p14 "covariance matrix
% 12.	 STATE p22 p23 p24
% 13.	 STATE p33 p34
% 14.	 STATE p44
% 15.	 NEW nysp1 ny1 nu1 nv1
% 16.	 NEW nth1 nth2 nth3 nth4
% 17.	 NEW nf1 nf2 nf3 nf4
% 18.	 NEW n11 n12 n13 n14 n22 n23 n24 n33 n34 n44
% 19.	 TIME t
% 20.	 TSAMP ts
% 21.	 INITIAL
% 22.	 "Compute sampled Am and Ao
% 23.	 a=exp(-z*w*h)
% 24.	 am1=-2*a*cos(w*h*sqrt(1-z*z))
% 25.
% 26.	 am2=a*a
% 27.	 aop=IF w*To>100 THEN 0 ELSE -exp(-h/To)
% 28.	 ao=IF cancel>0.5 THEN 0 ELSE -aop
% 29.	 SORT
% 30.	 "1.0 Parameter Estimation
% 31.	 "1.1 Computation of P*f and estimator gain k
% 32.	 pf1=p11*f1+p12*f2+p13*f3+p14*f4
% 33.	 pf2=p12*f1+p22*f2+p23*f3+p24*f4
% 34.	 pf3=p13*f1+p23*f2+p33*f3+p34*f4
% 35.	 pf4=p14*f1+p24*f2+p34*f3+p44*f4
% 36.	 denom=lambda+f1*pf1+f2*pf2+f3*pf3+f4*pf4
% 37.	 k1=pf1/denom
% 38.	 k2=pf2/denom
% 39.	 k3=pf3/denom
% 40.	 k4=pf4/denom
% 41.	 "1.2 Update estimates and covariances
% 42.	 eps=y-f1*th1-f2*th2-f3*th3-f4*th4
% 43.	 nth1=th1+k1*eps
% 44.	 nth2=th2+k2*eps
% 45.	 nth3=th3+k3*eps
% 46.	 nth4=th4+k4*eps
% 47.	 n11=(p11-pf1*k1)/lambda
% 48.	 n12=(p12-pf1*k2)/lambda
% 49.	 n13=(p13-pf1*k3)/lambda
% 50.	 n14=(p14-pf1*k4)/lambda
% 51.	 n22=(p22-pf2*k2)/lambda
% 52.	 n23=(p23-pf2*k3)/lambda
% 53.	 n24=(p24-pf2*k4)/lambda
% 54.	 n33=(p33-pf3*k3)/lambda
% 55.	 n34=(p34-pf3*k4)/lambda
% 56.	 n44=(p44-pf4*k4)/lambda
% 57.
% 58.	 "1.3 Update and filter regression vector
% 59.	 nf1=-y
% 60.	 nf2=f1
% 61.	 nf3=u
% 62.	 nf4=f3
% 63.	 "2.0 Control desig
% 64.	 b0=nth3
% 65.	 b1=nth4
% 66.	 "2.2 Solve the polynomial identity AR+BS=AoAm
% 67.	 n=b1*b1-a1*b0*b1+a2*b0*b0
% 68.	 r10=(ao*am2*b0^2+(a2-am2-ao*am1)*b0*b1+(ao+am1-a1)*b1^2)/n
% 69.	 w1=(a2*am1+a2*ao-a1*a2-am2*ao)*b0
% 70.	 s00=(w1+(-a1*am1-a1*ao-a2+a1^2+am2+am1*ao)*b1)/n
% 71.	 w2=(-a1*am2*ao+a2*am2+a2*am1*ao-a2^2)*b0
% 72.	 s10=(w2+(-a2*am1-a2*ao+a1*a2+am2*ao)*b1)/n
% 73.	 "2.3 Compute polynomial T=Ao*Am(1)/B(1)
% 74.	 bs=b0+b1
% 75.	 as=1+am1+am2
% 76.	 bm0=as/bs
% 77.	 "2.4 Choose control algorithm
% 78.	 r1=IF cancel>0.5 THEN b1/b0 ELSE r10
% 79.	 s0=IF cancel>0.5 THEN (am1-a1)/b0 ELSE s00
% 80.	 s1=IF cancel>0.5 THEN (am2-a2)/b0 ELSE s10
% 81.	 t0=IF cancel>0.5 THEN as/b0 ELSE bm0
% 82.	 t1=IF cancel>0.5 THEN 0 ELSE bm0*ao
% 83.	 "3.0 Control law with anti-windup
% 84.	 v=-ao*v1+t0*ysp+t1*ysp1-s0*y-s1*y1+(ao-r1)*u1
% 85.	 u=IF v<-ulim THEN -ulim ELSE IF v<ulim THEN v ELSE ulim
% 86.	 "3.1 Update controller state
% 87.	 ny1=y
% 88.	 nu1=u
% 89.	 nv1=v
% 90.	 nysp1=ysp
% 91.	 "4.0 Update sampling time
% 92.	 ts=t+h
% 93.	 "Parameters
% 94.	 lambda:1 "forgetting factor
% 95.	 To:200 "observer time constant
% 96.	 z:0.7 "desired closed loop damping
% 97.	 w:1 "desired closed loop natural frequency
% 98.	 h:1 "sampling period
% 99.	 ulim:1 "limit of control signal
% 100. cancel:1 "switch for cancellation
% 101. th1:-2 "initial estimates
% 102. 
% 103. th2:1
% 104. th3:0.01
% 105. th4:0.01
% 106. p11:100 "initial covariances
% 107. p22:100
% 108. p33:100
% 109. p44:100
% 110. END

%% Question 3.12
% The code for simulating Example 3.6 is listed below. Study the code and
% try to understand all the details. 
%%

% 1.	CONTINUOUS SYSTEM reg
% 2.	"Continuous time STR for the system b/[s(s+a)]
% 3.	"Desired response given by am2/(s^2+am1*s+am2)
% 4.	"Observer polynomial s+ao
% 5.	INPUT y ysp
% 6.	OUTPUT u
% 7.	STATE yf yf1 uf uf1 xu
% 8.	STATE th1 th2
% 9.	STATE p11 p12 p22
% 10.	DER dyf dyf1 duf duf1 dxu
% 11.	DER dth1 dth2
% 12.	DER dp11 dp12 dp22
% 13.	"Filter input and output
% 14.	dyf=yf1
% 15.	dyf1=-am1*yf1+am2*(y-yf)
% 16.	duf=uf1
% 17.	duf1=-am1*uf1+am2*(u-uf)
% 18.	"Update parameter estimate
% 19.	f1=-yf1
% 20.	f2=uf
% 21.	e=dyf1-f1*th1-f2*th2
% 22.	pf1=p11*f1+p12*f2
% 23.	pf2=p12*f1+p22*f2
% 24.	dth1=pf1*e
% 25.	dth2=pf2*e
% 26.	"Update covariance matrix
% 27.	dp11=alpha*p11-pf1*pf1
% 28.	dp12=alpha*p12-pf1*pf2
% 29.	dp22=alpha*p22-pf2*pf2
% 30.	det=p11*p22-p12*p12
% 31.	"Control design
% 32.	a=th1
% 33.	b=th2
% 34.	r1=am1+ao-a
% 35.	s0=(am2+am1*ao-a*r1)/b
% 36.	s1=am2*ao/b
% 37.	t0=am2/b
% 38.	"Control signal computation
% 39.	dxu=-ao*xu-(s1-ao*s0)*y+(ao-r1)*u
% 40.	v=t0*ysp-s0*y+xu
% 41.	u=if v<-ulim then -ulim else if v>ulim then ulim else v
% 42.	"Parameters
% 43.	am1:1.4
% 44.	am2:1
% 45.	alpha:0
% 46.	ao:2
% 47.	ulim:4
% 48.	END


%% Solution of the code of questions 3.5 and 3.12 
%%

% The descripition of the code is given as follows:
% Line 1 defines a DISCRETE TIME regulator. All the text after " is a comment.
% Inputs and outputs are defined with the reserved words INPUT and OUTPUT.
% All the variables between input and output are state variables with the
% reserved word STATE. From lines 7 to 14, all the variables are being
% defined. All the parameters that is being to be estimated are defined as
% with the type NEW. These parameters are defined from lines 15 to 20. The
% program execution takes place between INITIAL and END. Inside theses
% blocks, the variables are initialized and the estiamated parameters also.
% From lines 43 to 56, the covariance matrix is being updated recursively.
% Note the similarities withe these pseudo-code (?) with the implementaiton
% of the examples 3.4 and 3.5. Regarding the the question 3.12, the logic
% is the same of the example 3.6. Is important to note that the derivative
% filters are have the reserved word DER, but they are just  difference
% equations, once df(t)/dt ~= (f(t+h)-f(t))/h. 
