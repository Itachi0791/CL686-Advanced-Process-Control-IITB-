function K = KE(l1,lc1,l2,lc2,m1,m2,in7)
%KE
%    K = KE(L1,LC1,L2,LC2,M1,M2,IN7)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    19-Mar-2025 13:25:03

q1 = in7(1,:);
q2 = in7(2,:);
q1_d = in7(5,:);
q2_d = in7(6,:);
x0_d = in7(7,:);
y0_d = in7(8,:);
t2 = cos(q1);
t3 = sin(q1);
t4 = q1+q2;
t5 = q1_d+q2_d;
K = (m2.*((-x0_d+lc2.*t5.*sin(t4)+l1.*q1_d.*t3).^2+(y0_d+lc2.*t5.*cos(t4)+l1.*q1_d.*t2).^2))./2.0+(m1.*((x0_d-lc1.*q1_d.*t3).^2+(y0_d+lc1.*q1_d.*t2).^2))./2.0+(l1.^2.*m1.*q1_d.^2)./2.4e+1+(l2.^2.*m2.*t5.^2)./2.4e+1;
end
