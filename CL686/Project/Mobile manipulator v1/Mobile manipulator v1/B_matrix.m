function B_mat = B_matrix(x,u,params)
%B_matrix
%    B_mat = B_matrix(G,L1,LC1,L2,LC2,M1,M2,IN8,IN9)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    19-Mar-2025 13:25:10

g=params.g;l1=params.l1;lc1=params.lc1;l2=params.l2;lc2=params.lc2;m1=params.m1;m2=params.m2;

q1 = x(1,:);
q2 = x(2,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = sin(q1);
t5 = sin(q2);
t6 = q1+q2;
t7 = l1.^2;
t8 = l2.^2;
t9 = lc1.^2;
t10 = lc2.^2;
t11 = m1.^2;
t12 = m2.^2;
t13 = m2.^3;
t14 = q2.*2.0;
t19 = 1.0./m1;
t20 = 1.0./m2;
t21 = -q2;
t15 = cos(t14);
t16 = sin(t14);
t17 = cos(t6);
t18 = sin(t6);
t22 = -t14;
t23 = q2+t6;
t26 = q1+t21;
t30 = m1.*m2.*t8.*2.0;
t33 = t8.*t11;
t34 = t8.*t12;
t35 = t8.*t13;
t36 = m1.*m2.*t10.*1.2e+1;
t39 = t10.*t11.*1.2e+1;
t40 = m1.*t10.*t12.*1.2e+1;
t43 = l1.*lc1.*m1.*m2.*t8.*2.4e+1;
t44 = l1.*m1.*m2.*t2.*t8.*1.2e+1;
t45 = lc1.*m1.*m2.*t2.*t8.*1.2e+1;
t46 = l1.*m1.*m2.*t4.*t8.*1.2e+1;
t47 = lc1.*m1.*m2.*t4.*t8.*1.2e+1;
t49 = m1.*m2.*t7.*t8.*1.4e+1;
t50 = m1.*m2.*t8.*t9.*1.2e+1;
t52 = l1.*lc1.*m1.*m2.*t10.*1.44e+2;
t55 = l1.*m1.*m2.*t2.*t10.*7.2e+1;
t56 = lc1.*m1.*m2.*t2.*t10.*7.2e+1;
t59 = l1.*m1.*m2.*t4.*t10.*7.2e+1;
t60 = lc1.*m1.*m2.*t4.*t10.*7.2e+1;
t65 = m1.*m2.*t7.*t10.*8.4e+1;
t66 = m1.*m2.*t9.*t10.*7.2e+1;
t68 = lc1.*t2.*t10.*t11.*1.44e+2;
t69 = lc1.*t4.*t10.*t11.*1.44e+2;
t24 = cos(t23);
t25 = sin(t23);
t27 = q1+t22;
t28 = cos(t26);
t31 = sin(t26);
t37 = m2.*t33;
t38 = m1.*t34.*2.0;
t41 = m2.*t39;
t42 = t7.*t33;
t48 = l1.*lc1.*t34.*2.4e+1;
t51 = -t43;
t53 = l1.*t2.*t34.*1.2e+1;
t54 = lc1.*t2.*t33.*1.2e+1;
t57 = l1.*t4.*t34.*1.2e+1;
t58 = lc1.*t4.*t33.*1.2e+1;
t61 = t7.*t34.*1.3e+1;
t62 = t7.*t39;
t63 = t9.*t34.*1.2e+1;
t67 = -t52;
t70 = t15.*t52;
t71 = m1.*m2.*t7.*t10.*t15.*7.2e+1;
t78 = m1.*m2.*t9.*t10.*t15.*-7.2e+1;
t29 = cos(t27);
t32 = sin(t27);
t64 = -t48;
t73 = l1.*m1.*m2.*t10.*t24.*7.2e+1;
t74 = lc1.*m1.*m2.*t10.*t24.*7.2e+1;
t75 = l1.*m1.*m2.*t10.*t25.*7.2e+1;
t76 = lc1.*m1.*m2.*t10.*t25.*7.2e+1;
t77 = -t71;
t79 = -t73;
t80 = -t75;
t81 = t42+t49+t50+t51+t61+t62+t63+t64+t65+t66+t67+t70+t77+t78;
t82 = 1.0./t81;
mt1 = [0.0,0.0,0.0,0.0,t19.*t82.*(t30+t33+t34+t36+t39).*1.2e+1,t19.*t20.*t82.*(t35+t37+t38+t40+t41+l1.*lc2.*m1.*t3.*t12.*1.2e+1+l1.*lc2.*m2.*t3.*t11.*1.2e+1-lc1.*lc2.*m1.*t3.*t12.*1.2e+1-lc1.*lc2.*m2.*t3.*t11.*1.2e+1).*-1.2e+1,t19.*t82.*(t46+t47+t57+t58+t59+t60+t69+t76+t80),-t19.*t82.*(t44+t45+t53+t54+t55+t56+t68+t74+t79),0.0,0.0,0.0,0.0,t19.*t82.*(t30+t33+t34+t36+t39+l1.*lc2.*t3.*t11.*1.2e+1-lc1.*lc2.*t3.*t11.*1.2e+1+l1.*lc2.*m1.*m2.*t3.*1.2e+1-lc1.*lc2.*m1.*m2.*t3.*1.2e+1).*-1.2e+1];
mt2 = t19.*t20.*t82.*(t35+t37+t38+t40+t41+t7.*1.0./t19.^3+m1.*t7.*t12.*1.3e+1+m2.*t7.*t11.*1.4e+1+m1.*t9.*t12.*1.2e+1+m2.*t9.*t11.*1.2e+1-l1.*lc1.*m1.*t12.*2.4e+1-l1.*lc1.*m2.*t11.*2.4e+1+l1.*lc2.*m1.*t3.*t12.*2.4e+1+l1.*lc2.*m2.*t3.*t11.*2.4e+1-lc1.*lc2.*m1.*t3.*t12.*2.4e+1-lc1.*lc2.*m2.*t3.*t11.*2.4e+1).*1.2e+1;
mt3 = -t19.*t82.*(t46+t47+t57+t58+t59+t60+t69+t76+t80-lc2.*t7.*t11.*t18.*1.2e+1-lc2.*t9.*t11.*t18.*7.2e+1-lc2.*t9.*t11.*t31.*7.2e+1+l1.*lc1.*lc2.*t11.*t18.*7.2e+1+l1.*lc1.*lc2.*t11.*t31.*7.2e+1-lc2.*m1.*m2.*t7.*t18.*8.4e+1-lc2.*m1.*m2.*t9.*t18.*1.44e+2+lc2.*m1.*m2.*t7.*t31.*7.2e+1+l1.*lc1.*lc2.*m1.*m2.*t18.*2.16e+2-l1.*lc1.*lc2.*m1.*m2.*t31.*7.2e+1);
mt4 = [t19.*t82.*(t44+t45+t53+t54+t55+t56+t68+t74+t79-lc2.*t7.*t11.*t17.*1.2e+1-lc2.*t9.*t11.*t17.*7.2e+1-lc2.*t9.*t11.*t28.*7.2e+1+l1.*lc1.*lc2.*t11.*t17.*7.2e+1+l1.*lc1.*lc2.*t11.*t28.*7.2e+1-lc2.*m1.*m2.*t7.*t17.*8.4e+1-lc2.*m1.*m2.*t9.*t17.*1.44e+2+lc2.*m1.*m2.*t7.*t28.*7.2e+1+l1.*lc1.*lc2.*m1.*m2.*t17.*2.16e+2-l1.*lc1.*lc2.*m1.*m2.*t28.*7.2e+1),0.0,0.0,0.0,0.0,t19.*t82.*(l1.*m1.*m2.*t10.*t16.*6.0-lc1.*m1.*m2.*t10.*t16.*6.0).*-1.2e+1,t19.*t20.*t82.*(l1.*m1.*t10.*t12.*t16.*6.0+lc2.*m1.*t5.*t7.*t12.*1.3e+1+lc2.*m2.*t5.*t7.*t11+lc2.*m1.*t5.*t9.*t12.*1.2e+1-lc1.*m1.*t10.*t12.*t16.*6.0-l1.*lc1.*lc2.*m1.*t5.*t12.*2.4e+1).*1.2e+1];
mt5 = t19.*t82.*(t2.*t42+t2.*t50+t2.*t62+t2.*t66+m1.*m2.*t2.*t7.*t8.*1.3e+1+m1.*m2.*t2.*t7.*t10.*7.8e+1-m1.*m2.*t7.*t10.*t24.*4.2e+1-m1.*m2.*t9.*t10.*t24.*7.2e+1-m1.*m2.*t7.*t10.*t29.*3.6e+1-l1.*lc1.*m1.*m2.*t2.*t8.*2.4e+1-l1.*lc1.*m1.*m2.*t2.*t10.*1.44e+2+l1.*lc1.*m1.*m2.*t10.*t24.*1.08e+2+l1.*lc1.*m1.*m2.*t10.*t29.*3.6e+1);
mt6 = t19.*t82.*(t4.*t42+t4.*t50+t4.*t62+t4.*t66+m1.*m2.*t4.*t7.*t8.*1.3e+1+m1.*m2.*t4.*t7.*t10.*7.8e+1-m1.*m2.*t7.*t10.*t25.*4.2e+1-m1.*m2.*t9.*t10.*t25.*7.2e+1-m1.*m2.*t7.*t10.*t32.*3.6e+1-l1.*lc1.*m1.*m2.*t4.*t8.*2.4e+1-l1.*lc1.*m1.*m2.*t4.*t10.*1.44e+2+l1.*lc1.*m1.*m2.*t10.*t25.*1.08e+2+l1.*lc1.*m1.*m2.*t10.*t32.*3.6e+1);
B_mat = reshape([mt1,mt2,mt3,mt4,mt5,mt6],8,3);
end
