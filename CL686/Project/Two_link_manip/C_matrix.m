function C_matrix = C_matrix(in1,in2)
%C_matrix
%    C_matrix = C_matrix(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    14-Apr-2025 14:37:19

theta1 = in1(1,:);
theta2 = in1(2,:);
C_matrix = reshape([-sin(theta1),cos(theta1),sin(theta2).*(-4.0./5.0),cos(theta2).*(4.0./5.0),0.0,0.0,0.0,0.0],[2,4]);
end
