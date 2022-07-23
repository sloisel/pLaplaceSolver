function [ x,SOL ] = NesterovPathFollowingNew( c,F,x,varargin )
%[x,SOL] = NesterovPathFollowing( c,F,x )
%   The path-following algorithm described in the following paper:
%
%     Loisel, Sébastien. "Efficient algorithms for solving the p-Laplacian 
%     in polynomial time." Numerische Mathematik 146.2 (2020): 369-400.
%
%   This algorithm is an adaptive variant of Nesterov's path-following 
%   method to minimize c'*x subject to F<inf, wherein the size of the
%   t-steps is selected adaptively so that each centering iteration
%   converges in a reasonable number of damped Newton steps.
%
%   F should be a nu-self-concordant barrier on its domain and x
%   is the initial guess, and x should be feasible.
%
%   Optional parameters can be provided as follows:
%     'callback'  a callback function to be called at each iteration.
%                 callback function will be called as follows:
%                 callback(k,x,F0,G,H,t,PHASE);
%                 k is the iteration number; x is the iterate;
%                 F0,G,H are the function, gradient and Hessian; t is the
%                 barrier parameter and PHASE is either 'AUXILIARY' or
%                 'MAIN'.
%     'cholmode'  set to true to use chol, false for backslash. default is
%                 true.
%     'step'      growth factor for t barrier parameter. Defaults to 10. If
%                 you set 'step' to 1, you get the "short step" algorithm,
%                 of theoretical importance.
%     'maxit'     maximum number of iterations. Default is 30000.
%     'maxslow'   maximum number of "slow steps" before step-size
%                 adaptation is used. Defaults to 15. Set to inf to disable
%                 step size adaptation entirely.
%     'timeout'   If the iteration fails to converge before this number of
%                 seconds, the optimization fails with an error message.
%                 Default is 600 seconds.
%     'tol'       a tolerance. The outer iteration stops when t*tol>1.
%                 Default is 1e-6/length(c).
%     'verbose'   set to 1 to enable a callback function that prints all
%                 data.
%
%   Returns the converged iterate x, and SOL, a structure containing some
%   information about the iteration, including the values of t at each
%   iteration of the two phases of the algorithm.
SOL.maxit = 30000;
SOL.tol = 1e-6/length(c);
beta = 1/9;
gamma = 5/36;
SOL.stepsize = 10;
SOL.cholmode = true;
SOL.maxslow = 15;
SOL.checkderiv = false;
SOL.timeout = 600;
SOL.callback = @(k,x,F0,F1,F2,t,phase) 0;
SOL.F = F;
SOL.c = c;
SOL.x0 = x;
SOL.r = [];
    function N = starnorm(v,x)
        [~,~,H] = F(x);
        N = sqrt(v'*(H\v));
    end
SOL.starnorm = @starnorm;
for k=1:2:length(varargin)
    switch(varargin{k})
        case 'maxit'
            SOL.maxit = varargin{k+1};
        case 'tol'
            SOL.tol = varargin{k+1};
        case 'callback'
            SOL.callback = varargin{k+1};
        case 'verbose'
            SOL.callback = @myverbose;
        case 'step'
            SOL.stepsize = varargin{k+1};
        case 'cholmode'
            SOL.cholmode = varargin{k+1};
        case 'maxslow'
            SOL.maxslow = varargin{k+1};
        case 'timeout'
            SOL.timeout = varargin{k+1};
        case 'checkderiv'
            SOL.checkderiv = varargin{k+1};
        otherwise
            error('Unrecognized option');
    end
end
    function solver = makeSolver(H)
        if(strcmp(class(H),'function_handle'))
            internalSolver = H;
        else
            if(SOL.cholmode)
                [L,status,P] = chol(H,'lower','vector');
                if(status~=0)
                    warning('Turning chol off -- iteration will be slower');
                    SOL.cholmode = false;
                end
            end
            if(SOL.cholmode)
                Q(P) = 1:length(P);
                internalSolver = @cholSolver;
            else
                internalSolver = @(b) H\b;
            end
        end
        function y = cholSolver(b)
            y = L'\(L\b(P));
            y = y(Q);
        end
        function y = Hsol(b)
            assert(all(isreal(b)) && all(isfinite(b)),'b is not a real vector in H\b.');
            y = internalSolver(b);
            assert(all(isreal(y)) && all(isfinite(y)),'H\b is not a real vector. Probably H is numerically singular.');
        end
        solver = @Hsol;
    end
SOL.t = [1];
SOL.kappa = [SOL.stepsize];
SOL.accepted = [true];
slowdown = 0;
time0 = clock;
SOL.shortstep = (SOL.stepsize==1);
SOL.fastiter = 2;
for k=1:SOL.maxit
    assert(SOL.accepted(k) || ~SOL.shortstep(k),'Internal error: short steps should always be accepted.');
    [F0,G,H] = F(x);
    SOL.callback(k,x,F0,G,H,SOL.t(k),'AUXILIARY');
    assert(isfinite(F0) && isreal(F0),'F(x) is not a real number');
    if(SOL.checkderiv)
        [G1,H1] = mydiff(F,x);
        eG = norm(G1-G);
        tG = SOL.checkderiv*(abs(F0)+norm(G1));
        assert(~(eG>tG),'Error in gradient of %e exceeds tolerance %e',eG,tG);
        eH = norm(H1-H,inf)/norm(H1,inf);
        tH = SOL.checkderiv*(abs(F0)+norm(G1)+norm(H1,inf));
        assert(~(eH>tH),'Error in Hessian of %e exceeds tolerance %e',eH,tH);
    end
    if(k==1), G0 = G; end
    sol = makeSolver(H);
    HsolG = sol(G);
    if(sqrt(G'*HsolG)<=sqrt(beta)/(1+sqrt(beta)))
        x = x - HsolG;
        break;
    end
    HsolG0 = sol(G0);
    NN = sqrt(dot(-SOL.t(k)*G0+G,-SOL.t(k)*HsolG0+HsolG));
    SOL.lambda(k) = NN;
    SOL.accepted(k+1) = false;
    SOL.kappa(k+1) = SOL.kappa(k);
    SOL.shortstep(k+1) = false;
    SOL.t(k+1) = SOL.t(k);
    if(NN<0.2 || SOL.accepted(k))
        SOL.accepted(k) = true;
        if(slowdown>SOL.maxslow/2)
            SOL.kappa(k+1) = sqrt(SOL.kappa(k));
        elseif(slowdown<=SOL.fastiter && k>1)
            SOL.kappa(k+1) = min(SOL.stepsize,SOL.kappa(k)^2);
        end
        slowx = x;
        slowt = SOL.t(k);
        slowdown = 0;
        medt = SOL.t(k)-gamma/sqrt(G0'*HsolG0);
        fastt = SOL.t(k)/SOL.kappa(k+1);
        shortstep = (medt<=fastt);
        SOL.t(k+1) = min(medt,fastt);
        if(shortstep)
            SOL.shortstep(k) = true;
            SOL.accepted(k+1) = true;
        end
    else
        slowdown = slowdown+1;
        if(slowdown == SOL.maxslow)
            SOL.accepted(k+1) = true;
            SOL.t(k+1) = slowt;
            x = slowx;
            SOL.kappa(k+1) = (SOL.kappa(k))^0.25;
            slowdown = SOL.fastiter+1;
            continue;
        end
    end
    t = SOL.t(k+1);
    del = t*HsolG0-HsolG;
    if(SOL.shortstep(k))
        r = 1;
    else
        t = SOL.t(k+1);
        r = backtracking(@(r) -t*dot(G0,x+r*del)+F(x+r*del),...
                         -t*dot(G0,x)+F0, ...
                         0.01*dot(-t*G0+G,del));
    end
    SOL.r(k) = r;
    x = x + r*del;
    time1 = clock;
    SOL.elapsed = etime(time1,time0);
    assert(SOL.elapsed<=SOL.timeout,'Exceeded maximum running time in auxiliary phase');
end
assert(k<SOL.maxit,'Iteration count exceeded during auxiliary path following algorithm');
k0 = k+1;
SOL.t(k0) = 0;
SOL.kappa(k0) = SOL.kappa(end);
SOL.accepted(k0) = true;
SOL.shortstep(k0) = (SOL.stepsize==1);
slowdown = 0;
for k=k0:SOL.maxit
    assert(SOL.accepted(k) || ~SOL.shortstep(k),'Internal error: short steps should always be accepted.');
    [F0,G,H] = F(x);
    SOL.callback(k,x,F0,G,H,SOL.t(k),'MAIN');
    assert(isfinite(F0) && isreal(F0),'F(x) is not a real number');
    if(SOL.checkderiv)
        [G1,H1] = mydiff(F,x);
        eG = norm(G1-G);
        tG = SOL.checkderiv*(abs(F0)+norm(G1));
        eH = norm(H1-H,inf)/norm(H1,inf);
        tH = SOL.checkderiv*(abs(F0)+norm(G1)+norm(H1,inf));
    end
%    Hsolc = Hsol(H,c);
%    HsolG = Hsol(H,G);
    sol = makeSolver(H);
    Hsolc = sol(c);
    HsolG = sol(G);
    t = SOL.t(k);
    del = t*Hsolc+HsolG;
    NN = sqrt((t*c+G)'*del);
    SOL.lambda(k) = NN;
    SOL.accepted(k+1) = false;
    SOL.kappa(k+1) = SOL.kappa(k);
    SOL.shortstep(k+1) = false;
    SOL.t(k+1) = SOL.t(k);
    if(NN<0.2 || SOL.accepted(k))
        SOL.accepted(k) = true;
        if(SOL.tol*t>=1)
            break;
        end
        if(slowdown>SOL.maxslow/2)
            SOL.kappa(k+1) = sqrt(SOL.kappa(k));
        elseif(SOL.fastiter<=1 && k>k0+1)
            SOL.kappa(k+1) = min(SOL.stepsize,SOL.kappa(k)^2);
        end
        slowx = x;
        slowt = t;
        slowdown = 0;
        medt = SOL.t(k) + gamma/sqrt(c'*Hsolc);
        fastt = SOL.t(k)*SOL.kappa(k+1);
        shortstep = medt>=fastt;
        SOL.t(k+1) = max(medt,fastt);
        if(shortstep)
            SOL.shortstep(k) = true;
            SOL.accepted(k+1) = true;
        end
    else
        slowdown = slowdown+1;
        if(slowdown == SOL.maxslow)
            SOL.accepted(k+1) = true;
            SOL.t(k+1) = slowt;
            x = slowx;
            SOL.kappa(k+1) = (SOL.kappa(k))^0.25;
            slowdown = SOL.fastiter+1;
            continue;
        end
    end
    t = SOL.t(k+1);
    del = t*Hsolc+HsolG;
    if(SOL.shortstep(k))
        r = 1;
    else
        r = backtracking(@(r) t*dot(c,x-r*del)+F(x-r*del), ...
                         t*dot(c,x)+F0, ...
                         0.01*dot(t*c+G,-del));
    end
    SOL.r(k) = r;
    x = x-r*del;
    time1 = clock;
    SOL.elapsed = etime(time1,time0);
    assert(SOL.elapsed<=SOL.timeout,'Exceeded maximum running time in auxiliary phase');
end
assert(k<SOL.maxit,'Iteration count exceeded during main path following algorithm');
SOL.x = x;
end

function t = backtracking(F,F0,a)
t = 1;
F1 = inf;
while(true)
    Fprev = F1;
    F1 = F(t);
    if(isreal(F1) && isfinite(F1) && (F1<=F0+a*t || F1>=Fprev))
        return;
    end
    t = 0.25*t;
    assert(t>0,'Backtracking search error: the Newton direction is not a descent direction');
end
end

function myverbose(k,x,F0,F1,F2,t,phase)
disp(sprintf('%s: Iteration %d',phase,k));
x
t
F0
F1
F2
end

