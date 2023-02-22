% test 6pt


function  [F,f] = solve_f_E_f_elimIdeal(x,xp,T1,T2,met)

F{1}= 10000*eye(3,3);
f(1) = 10000;

if nargin<3
    T1 = eye(3);
    T2 = eye(3);
	met = 2;
end

if nargin<5
    met = 2;
end
    
%fill the 6x9 matrix
% cplumns correspond to monomials
% f11, f12, f13, f21, f22, f23, f31, f32,  f33
M = [ x(:,1).*xp(:,1), x(:,1).*xp(:,2), x(:,1).*xp(:,3), x(:,2).*xp(:,1), x(:,2).*xp(:,2), x(:,2).*xp(:,3), x(:,3).*xp(:,1), x(:,3).*xp(:,2), x(:,3).*xp(:,3)];

%null space
n = null(M);

%GB solver
[a, b] = solver_sw6pt_sturmfels(n(:,1), n(:,2), n(:,3));

% reconstruct solutions
for i = 1:length(a)
    % 9dim vector
    FF = a(i)*n(:,1)+b(i)*n(:,2)+n(:,3);
    
    %fundamental matrix
    F{i} = T1'*reshape(FF(1:9)', 3,3)*T2 ;
    
    x11 = FF(1);
    x12 = FF(2);
    x13 = FF(3);
    x21 = FF(4);
    x22 = FF(5);
    x23 = FF(6);
    x31 = FF(7);
    x32 = FF(8);
    x33 = FF(9);
    
    if (met == 1)
        [U D V] = svd(F{i});
        u1 = U(:,1);
        u2 = U(:,2);
        v1 = V(:,1);
        v2 = V(:,2);
        as = D(1,1);
        bs = D(2,2);
        %f(i) =  sqrt((-(u2(3)*v1(3)*(as*u1(3)*v1(3)+bs*u2(3)*v2(3)))/(as*u1(3)*u2(3)*(1-v1(3)^2)+bs*v1(3)*v2(3)*(1-u2(3)^2))));
        f(i) = sqrt((-x13^2*x32*x33-x23^2*x32*x33+x12*x13*x33^2+x22*x23*x33^2) /(x11*x13*x31*x32+x21*x23*x31*x32+x12*x13*x32^2+x22*x23*x32^2-x11*x12*x31*x33-x21*x22*x31*x33-x12^2*x32*x33-x22^2*x32*x33));
    end
    
    if (met == 2)
        f(i) = Solve_f(F{i});
    end
    
    if (met == 3)
        f2 = F2f1f2(F{i});
        f(i) = 1/2*(f2(1)+f2(2));
    end
    
    if (met == 4)
    	f(i) = f_fFf(F{i});
    end

end
 
