function B_mat = B_matrix(x,u,params)
%B_matrix
%    B_mat = B_matrix(G,L1,LC1,L2,LC2,M1,M2,IN8,IN9)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    20-Mar-2025 13:27:31
g=params.g;l1=params.l1;lc1=params.lc1;l2=params.l2;lc2=params.lc2;m1=params.m1;m2=params.m2;

q1 = x(1,:);
q2 = x(2,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = sin(q1);
t5 = q1+q2;
t6 = l1.^2;
t7 = l2.^2;
t8 = lc1.^2;
t9 = lc2.^2;
t10 = m1.^2;
t11 = m2.^2;
t12 = m2.^3;
t13 = q1.*2.0;
t14 = q2.*2.0;
t20 = 1.0./m1;
t21 = 1.0./m2;
t22 = -q2;
t15 = cos(t13);
t16 = cos(t14);
t17 = sin(t13);
t18 = cos(t5);
t19 = sin(t5);
t23 = q2+t5;
t26 = q1+t22;
t27 = t5.*2.0;
t29 = m1.*m2.*t7.*2.0;
t31 = t7.*t10;
t32 = t7.*t11;
t33 = t7.*t12;
t34 = m1.*m2.*t9.*1.2e+1;
t39 = t9.*t10.*1.2e+1;
t40 = m1.*t9.*t11.*1.2e+1;
t43 = l1.*lc1.*m1.*m2.*t7.*1.2e+1;
t44 = l1.*lc1.*m1.*m2.*t7.*2.4e+1;
t45 = l1.*m1.*m2.*t2.*t7.*1.2e+1;
t46 = lc1.*m1.*m2.*t2.*t7.*1.2e+1;
t47 = l1.*m1.*m2.*t4.*t7.*1.2e+1;
t48 = lc1.*m1.*m2.*t4.*t7.*1.2e+1;
t52 = m1.*m2.*t6.*t7.*1.3e+1;
t53 = m1.*m2.*t6.*t7.*1.4e+1;
t54 = m1.*m2.*t7.*t8.*1.2e+1;
t57 = l1.*lc1.*m1.*m2.*t9.*7.2e+1;
t58 = l1.*lc1.*m1.*m2.*t9.*1.44e+2;
t61 = l1.*m1.*m2.*t2.*t9.*7.2e+1;
t62 = lc1.*m1.*m2.*t2.*t9.*7.2e+1;
t65 = l1.*m1.*m2.*t4.*t9.*7.2e+1;
t66 = lc1.*m1.*m2.*t4.*t9.*7.2e+1;
t71 = m1.*m2.*t6.*t9.*7.8e+1;
t72 = m1.*m2.*t6.*t9.*8.4e+1;
t73 = m1.*m2.*t8.*t9.*7.2e+1;
t76 = lc1.*t2.*t9.*t10.*1.44e+2;
t77 = lc1.*t4.*t9.*t10.*1.44e+2;
t78 = t8.*t9.*t10.*7.2e+1;
t24 = cos(t23);
t25 = sin(t23);
t28 = cos(t26);
t30 = sin(t26);
t35 = cos(t27);
t36 = m2.*t31;
t37 = sin(t27);
t38 = m1.*t32.*2.0;
t41 = m2.*t39;
t42 = t6.*t31;
t49 = t6.*t32.*6.0;
t50 = t8.*t31.*6.0;
t51 = l1.*lc1.*t32.*2.4e+1;
t55 = -t43;
t56 = -t44;
t59 = l1.*t2.*t32.*1.2e+1;
t60 = lc1.*t2.*t31.*1.2e+1;
t63 = l1.*t4.*t32.*1.2e+1;
t64 = lc1.*t4.*t31.*1.2e+1;
t67 = t6.*t32.*1.3e+1;
t68 = t6.*t39;
t69 = t8.*t32.*1.2e+1;
t74 = -t57;
t75 = -t58;
t79 = t15.*t43;
t80 = t17.*t43;
t83 = t15.*t57;
t84 = t16.*t57;
t85 = t16.*t58;
t86 = t17.*t57;
t89 = m1.*m2.*t6.*t9.*t16.*7.2e+1;
t96 = m1.*m2.*t8.*t9.*t16.*-7.2e+1;
t97 = t17.*t78;
t100 = t15.*t78;
t70 = -t51;
t81 = t17.*t49;
t82 = t17.*t50;
t87 = t15.*t49;
t88 = t15.*t50;
t91 = l1.*m1.*m2.*t9.*t24.*7.2e+1;
t92 = lc1.*m1.*m2.*t9.*t24.*7.2e+1;
t93 = l1.*m1.*m2.*t9.*t25.*7.2e+1;
t94 = lc1.*m1.*m2.*t9.*t25.*7.2e+1;
t95 = -t89;
t101 = m1.*m2.*t6.*t9.*t35.*6.0;
t102 = m1.*m2.*t6.*t9.*t37.*6.0;
t103 = t35.*t57;
t105 = t35.*t73;
t106 = t37.*t73;
t107 = l1.*lc1.*m1.*m2.*t9.*t37.*-7.2e+1;
t98 = -t91;
t99 = -t93;
t108 = t80+t81+t82+t86+t97+t102+t106+t107;
t109 = t42+t53+t54+t56+t67+t68+t69+t70+t72+t73+t75+t85+t95+t96;
t110 = 1.0./t109;
t111 = t20.*t108.*t110;
t112 = -t111;
mt1 = [0.0,0.0,0.0,0.0,t20.*t110.*(t29+t31+t32+t34+t39).*1.2e+1,t20.*t21.*t110.*(t33+t36+t38+t40+t41+l1.*lc2.*m1.*t3.*t11.*1.2e+1+l1.*lc2.*m2.*t3.*t10.*1.2e+1-lc1.*lc2.*m1.*t3.*t11.*1.2e+1-lc1.*lc2.*m2.*t3.*t10.*1.2e+1).*-1.2e+1,t20.*t110.*(t47+t48+t63+t64+t65+t66+t77+t94+t99),-t20.*t110.*(t45+t46+t59+t60+t61+t62+t76+t92+t98),0.0,0.0,0.0,0.0,t20.*t110.*(t29+t31+t32+t34+t39+l1.*lc2.*t3.*t10.*1.2e+1-lc1.*lc2.*t3.*t10.*1.2e+1+l1.*lc2.*m1.*m2.*t3.*1.2e+1-lc1.*lc2.*m1.*m2.*t3.*1.2e+1).*-1.2e+1];
mt2 = [t20.*t21.*t110.*(t33+t36+t38+t40+t41+t6.*1.0./t20.^3+m1.*t6.*t11.*1.3e+1+m2.*t6.*t10.*1.4e+1+m1.*t8.*t11.*1.2e+1+m2.*t8.*t10.*1.2e+1-l1.*lc1.*m1.*t11.*2.4e+1-l1.*lc1.*m2.*t10.*2.4e+1+l1.*lc2.*m1.*t3.*t11.*2.4e+1+l1.*lc2.*m2.*t3.*t10.*2.4e+1-lc1.*lc2.*m1.*t3.*t11.*2.4e+1-lc1.*lc2.*m2.*t3.*t10.*2.4e+1).*1.2e+1];
mt3 = [-t20.*t110.*(t47+t48+t63+t64+t65+t66+t77+t94+t99-lc2.*t6.*t10.*t19.*1.2e+1-lc2.*t8.*t10.*t19.*7.2e+1-lc2.*t8.*t10.*t30.*7.2e+1+l1.*lc1.*lc2.*t10.*t19.*7.2e+1+l1.*lc1.*lc2.*t10.*t30.*7.2e+1-lc2.*m1.*m2.*t6.*t19.*8.4e+1-lc2.*m1.*m2.*t8.*t19.*1.44e+2+lc2.*m1.*m2.*t6.*t30.*7.2e+1+l1.*lc1.*lc2.*m1.*m2.*t19.*2.16e+2-l1.*lc1.*lc2.*m1.*m2.*t30.*7.2e+1)];
mt4 = [t20.*t110.*(t45+t46+t59+t60+t61+t62+t76+t92+t98-lc2.*t6.*t10.*t18.*1.2e+1-lc2.*t8.*t10.*t18.*7.2e+1-lc2.*t8.*t10.*t28.*7.2e+1+l1.*lc1.*lc2.*t10.*t18.*7.2e+1+l1.*lc1.*lc2.*t10.*t28.*7.2e+1-lc2.*m1.*m2.*t6.*t18.*8.4e+1-lc2.*m1.*m2.*t8.*t18.*1.44e+2+lc2.*m1.*m2.*t6.*t28.*7.2e+1+l1.*lc1.*lc2.*m1.*m2.*t18.*2.16e+2-l1.*lc1.*lc2.*m1.*m2.*t28.*7.2e+1),0.0,0.0,0.0,0.0,t20.*t110.*(l1.*t4.*t32+lc1.*t4.*t31+lc1.*t4.*t39+l1.*m1.*m2.*t4.*t7+l1.*m1.*m2.*t4.*t9.*6.0-l1.*m1.*m2.*t9.*t25.*6.0+lc1.*m1.*m2.*t4.*t7+lc1.*m1.*m2.*t4.*t9.*6.0+lc1.*m1.*m2.*t9.*t25.*6.0).*1.2e+1];
mt5 = [t20.*t21.*t110.*(l1.*t4.*t33+lc1.*t4.*t36+lc1.*t4.*t41+l1.*m1.*t4.*t32+lc1.*m1.*t4.*t32+l1.*m1.*t4.*t9.*t11.*6.0-l1.*m1.*t9.*t11.*t25.*6.0+lc1.*m1.*t4.*t9.*t11.*6.0-lc2.*m1.*t6.*t11.*t19.*7.0-lc2.*m2.*t6.*t10.*t19-lc2.*m1.*t8.*t11.*t19.*1.2e+1-lc2.*m2.*t8.*t10.*t19.*6.0+lc1.*m1.*t9.*t11.*t25.*6.0+lc2.*m1.*t6.*t11.*t30.*6.0-lc2.*m2.*t8.*t10.*t30.*6.0+l1.*lc1.*lc2.*m1.*t11.*t19.*1.8e+1+l1.*lc1.*lc2.*m2.*t10.*t19.*6.0-l1.*lc1.*lc2.*m1.*t11.*t30.*6.0+l1.*lc1.*lc2.*m2.*t10.*t30.*6.0).*-1.2e+1];
mt6 = [t20.*t110.*(t42+t49+t50+t52+t54+t55+t68+t71+t73+t74+t78+t84+t95-t101+t103-t6.*t15.*t32.*6.0-t8.*t15.*t31.*6.0-t8.*t9.*t10.*t15.*7.2e+1-m1.*m2.*t8.*t9.*t35.*7.2e+1-l1.*lc1.*m1.*m2.*t7.*t15.*1.2e+1-l1.*lc1.*m1.*m2.*t9.*t15.*7.2e+1),t112,0.0,0.0,0.0,0.0,t20.*t110.*(l1.*t2.*t32+lc1.*t2.*t31+lc1.*t2.*t39+l1.*m1.*m2.*t2.*t7+l1.*m1.*m2.*t2.*t9.*6.0-l1.*m1.*m2.*t9.*t24.*6.0+lc1.*m1.*m2.*t2.*t7+lc1.*m1.*m2.*t2.*t9.*6.0+lc1.*m1.*m2.*t9.*t24.*6.0).*-1.2e+1];
mt7 = [t20.*t21.*t110.*(l1.*t2.*t33+lc1.*t2.*t36+lc1.*t2.*t41+l1.*m1.*t2.*t32+lc1.*m1.*t2.*t32+l1.*m1.*t2.*t9.*t11.*6.0-l1.*m1.*t9.*t11.*t24.*6.0+lc1.*m1.*t2.*t9.*t11.*6.0-lc2.*m1.*t6.*t11.*t18.*7.0-lc2.*m2.*t6.*t10.*t18-lc2.*m1.*t8.*t11.*t18.*1.2e+1-lc2.*m2.*t8.*t10.*t18.*6.0+lc1.*m1.*t9.*t11.*t24.*6.0+lc2.*m1.*t6.*t11.*t28.*6.0-lc2.*m2.*t8.*t10.*t28.*6.0+l1.*lc1.*lc2.*m1.*t11.*t18.*1.8e+1+l1.*lc1.*lc2.*m2.*t10.*t18.*6.0-l1.*lc1.*lc2.*m1.*t11.*t28.*6.0+l1.*lc1.*lc2.*m2.*t10.*t28.*6.0).*1.2e+1,t112,t20.*t110.*(t42+t49+t50+t52+t54+t55+t68+t71+t73+t74+t78+t79+t83+t84+t87+t88+t95+t100+t101+t105-l1.*lc1.*m1.*m2.*t9.*t35.*7.2e+1)];
B_mat = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7],8,4);
end
