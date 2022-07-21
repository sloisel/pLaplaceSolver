function F = pLaplaceBarrierY( p,D,DG,g,R,omega)
d = length(D);
m = size(D{1},1);
n = size(D{1},2);
pinf = false;
if(p==inf)
    pinf = true;
    p = 1;
end
if(p<2), sigma = 2;
else sigma = 1;
end
cholmode = true;
OMEGA = sum(omega);
    function [F0,F1,F2] = barrier(x)
        u = x(1:n);
        s = x(n+1:end);
        if(pinf)
            s = s*ones(m,1);
        end
        y = cell(1,d);
        su = zeros(m,1);
        for j=1:d
            y{j} = D{j}*u+DG{j}*g;
            su = su+y{j}.^2;
        end
        z = s.^(2/p)-su;
        tau = R - omega.*s;
        F0 = -sum(logp(z))-sigma*sum(logp(s))-sum(logp(tau));
%        assert(isreal(F0) && isfinite(F0));
        if(nargout==1), return; end

        Fu = zeros(n,1);
        for j=1:d
            Fu = Fu + 2*(D{j})'*(y{j}./z);
        end
        Fs = -(2/p).*(1./z).*s.^((2/p)-1)-sigma./s+(omega./tau);
        F1 = [Fu;Fs];
        if(pinf)
            EE = blkdiag(speye(n),sparse(ones(1,m)));
            F1 = EE*F1;
        end
        if(nargout==2), return; end
        
        Z1 = spdiags(z.^-1,0,m,m);
        Z2 = spdiags(z.^-2,0,m,m);
        E = sparse(n,n);
        SU = sparse(m,n);
        for j=1:d
            YjDj = spdiags(y{j},0,m,m)*(D{j});
            SU = SU + spdiags(y{j},0,m,m)*(D{j});
            E = E + 2*(D{j})'*Z1*(D{j});
                  %+ 4*YjDj'*Z2*YjDj;
        end
        E = E + 4*SU'*Z2*SU;
        Sp = spdiags(s.^(2/p-1),0,m,m);
        K = (-4/p)*SU'*Z2*Sp;
        G = -(2/p)*(2/p-1)*Z1*spdiags(s.^(2/p-2),0,m,m) + ...
            (4/p^2)*Z2*spdiags(s.^(4/p-2),0,m,m) + ...
            sigma*spdiags(s.^-2,0,m,m) + ...
            spdiags((omega./tau).^2,0,m,m);
        H = [E K;K' G];
        if(pinf)
            H = EE*H*(EE');
        end
        F2 = H;
    end
F = @barrier;
end

function y = logp(x)
y = log(max(x,0));
end

