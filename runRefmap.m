function [model, llh] = runRefmap(Z, S, Ld, N, Ik, q)
%% Run VI algorithm for RefMap (REgional Fine-MAPping)
%  Author: Sai Zhang (zhangsai@stanford.edu)
%
% N: number of GWAS samples
% I: list of SNP number in each LD, K x 1
% J: list of region number in each LD, K x 1
% K: number of LD blocks
% Ik: cell of SNP number in each region of LD, {I1, ..., IK}, Ik: J_k x 1
% Z: cell of z-score vectors, {Z_1, ..., Z_K}, Z_k: I_k x 1
% S: cell of annotations, {A_1, ..., A_K}, A_k: J_k x D
% D: number of epigenetic features + bias term
% Ld: cell of LD matrix, {LD_1, ..., LD_K}, LD_k: I_k x I_k


K = length(Z);
I = [];
J = [];
for k=1:K
    I = [I, size(Z{k},1)];
    J = [J, size(S{k},1)];
end
D = size(S{1},2);

fprintf('Running VB-EM for RefMap...\n');
tol = 1e-6;
maxiter = 1000;
llh = -inf(1, maxiter);

model = init(Z, S, Ik, J, D, K, Ld, q);
for iter = 2:maxiter
    %E-step
    model = expect(Z, S, Ld, Ik, I, D, K, N, model);
    %M-step
    model = maximize(S, K, model);
    %eval
    llh(iter) = eval(S, Z, Ld, Ik, K, N, model);
    fprintf('Iter %d: llh is %6.4f\n', iter-1, llh(iter));
    
    if abs(llh(iter)-llh(iter-1)) < tol*abs(llh(iter))
        fprintf('Converged!\n')
        break; 
    end
end

llh = llh(2:iter);
end

function s = sigmoid(x)
%%
s = 1./(1 + exp(-x));
end

function f = fchi(x)
%%
f = (0.5./x) .* (sigmoid(x)-0.5);
end

function model = init(Z, S, Ik, J, D, K, Ld, q)
%% Model initialization

model.mu0 = 0;
model.beta0 = 1;
model.a00 = 400;
model.b00 = 1;
model.a01 = 1e-6;
model.b01 = 1e-6;
model.nu0 = D;
model.W0 = eye(D);
model.W0inv = eye(D);

model.lnLd = zeros(K, 1);
for k=1:K
    model.q = q;
    
    U = chol(Ld{k});
    model.lnLd(k) = 2*sum(log(diag(U)));
    
    model.Eqmk{k} = zeros([J(k),1]);
    model.Eqlambdak{k} = zeros([J(k),1]);
    for j=1:J(k)
        model.Eqmk{k}(j) = mean(Z{k}( sum(Ik{k}(1:j))-Ik{k}(j)+1 : sum(Ik{k}(1:j)) ));
        model.Eqlambdak{k}(j) = 1/var(Z{k});
    end

    y = zeros([J(k),3]);
    y(:,1) = gamrnd(model.q/2,1,[J(k),1]);
    y(:,2) = gamrnd(model.q/2,1,[J(k),1]);
    y(:,3) = gamrnd(1-model.q,1,[J(k),1]);
    model.Eqtkneg{k} = y(:,1)./sum(y,2).*(sum(S{k}(:,1:1),2)~=0);
    model.Eqtkpos{k} = y(:,2)./sum(y,2).*(sum(S{k}(:,1:1),2)~=0);
    model.Eqtk0{k} = 1-model.Eqtkneg{k}-model.Eqtkpos{k};
    
    model.EqmkVec{k} = [];
    model.EqLambdakVec{k} = [];
    for j=1:length(Ik{k})
        model.EqmkVec{k} = [model.EqmkVec{k}; model.Eqmk{k}(j)*ones([Ik{k}(j),1])];
        model.EqLambdakVec{k} = [model.EqLambdakVec{k}; model.Eqlambdak{k}(j)*ones([Ik{k}(j),1])];
    end
end

model.Eqtauneg = gamrnd(1,1);
model.Eqtau0 = gamrnd(1,1);
model.Eqtaupos = gamrnd(1,1);
model.Eqvneg = gamrnd(5,1);
model.Vqvneg = gamrnd(1,1);
model.Eqvneg2 = model.Vqvneg + model.Eqvneg^2;
model.Eqvpos = gamrnd(5,1);
model.Vqvpos = gamrnd(1,1);
model.Eqvpos2 = model.Vqvpos + model.Eqvpos^2;

model.Eqmneg = normrnd(0,1);
model.Eqmpos = normrnd(0,1);
model.Eqlambdaneg = gamrnd(1,1);
model.Eqlambdapos = gamrnd(1,1);
model.Eqw = normrnd(0,1,[D,1]);

model.Eqw(end) = -log(1/model.q-1) - sum(model.Eqw(1:end-1));
model.Vqw = wishrnd(eye(D),D);
model.Eqw2 = model.Vqw + model.Eqw*model.Eqw';
model.EqLambda = wishrnd(eye(D),D);

%init xi
model = maximize(S, K, model);
end

function model = expect(Z, S, Ld, Ik, I, D, K, N, model)
%% E-step
% 1. q(u_k)
for k=1:K
    model.tildeLambdak{k} = N*Ld{k} + diag(model.EqLambdakVec{k});
    U = chol(model.tildeLambdak{k});
    tmp = sqrt(N)*Z{k} + diag(model.EqLambdakVec{k})*model.EqmkVec{k};
    model.tildemuk{k} = U \ (U' \ tmp);
    
    model.Equk{k} = model.tildemuk{k}; %I_k x 1
    model.Vquk{k} = U \ (U' \ eye(I(k))); %I_k x I_k
    model.Equk2{k} = model.Vquk{k} + model.Equk{k}*model.Equk{k}'; %I_k x I_k
    %model.Hquk{k} = -sum(log(diag(U))) + 0.5*I(k)*(1+log(2*pi));
    
    for j=1:length(Ik{k})
        model.Equksum{k}(j,1) = sum(model.Equk{k}( (sum(Ik{k}(1:j))-Ik{k}(j)+1) : sum(Ik{k}(1:j)) )); %J_k x 1, sum_i{E[u_{ijk}]}
        tmp = diag(model.Equk2{k}); %I_k x 1
        model.Equk2sum{k}(j,1) = sum(tmp( (sum(Ik{k}(1:j))-Ik{k}(j)+1) : sum(Ik{k}(1:j)) )); %J_k x 1, sum_i{E[u_{ijk}^2]}
    end
end

% 2. q(m_{j,k})
for k=1:K
    model.tildelambdak{k} = Ik{k}.*model.Eqlambdak{k} + model.Eqtkneg{k}*model.Eqtauneg + model.Eqtk0{k}*model.Eqtau0 + model.Eqtkpos{k}*model.Eqtaupos;
    tmp = model.Eqlambdak{k}.*model.Equksum{k} - model.Eqvneg*model.Eqtauneg*model.Eqtkneg{k} + model.Eqvpos*model.Eqtaupos*model.Eqtkpos{k};
    model.tildemuk{k} = tmp .* model.tildelambdak{k}.^(-1);
    
    model.Eqmk{k} = model.tildemuk{k}; %J_k x 1
    model.Vqmk{k} = 1./model.tildelambdak{k}; %J_k x 1
    model.Eqmk2{k} = model.Vqmk{k} + model.Eqmk{k}.^2; %J_k x 1
    %model.Hqmk{k} = 0.5*log(model.Vqmk{k}) + 0.5*(1+log(2*pi));
    
    model.EqmkVec{k} = [];
    for j=1:length(Ik{k})
        model.EqmkVec{k} = [model.EqmkVec{k}; model.Eqmk{k}(j)*ones([Ik{k}(j),1])];
    end
end

% 3. q(lambda_{j,k})
for k=1:K
    model.tildeak{k} = model.a00 + 0.5*Ik{k};
    model.tildebk{k} = model.b00 + 0.5*model.Equk2sum{k} + 0.5*Ik{k}.*model.Eqmk2{k} - model.Eqmk{k}.*model.Equksum{k};
    
    model.Eqlambdak{k} = model.tildeak{k} ./ model.tildebk{k};
    model.Vqlambdak{k} = model.tildeak{k} ./ (model.tildebk{k}).^2;
    model.Eqlnlambdak{k} = psi(model.tildeak{k}) - log(model.tildebk{k});
    %model.Hqlambdak{k} = gammaln(model.tildeak{k}) - (model.tildeak{k}-1).*psi(model.tildeak{k}) - log(model.tildebk{k}) + model.tildeak{k};
    
    model.EqLambdakVec{k} = [];
    for j=1:length(Ik{k})
        model.EqLambdakVec{k} = [model.EqLambdakVec{k}; model.Eqlambdak{k}(j)*ones([Ik{k}(j),1])];
    end
end

% 4. q(tau_{-1})
tmp1 = 0;
tmp2 = 0;
for k=1:K
    tmp1 = tmp1 + sum(model.Eqtkneg{k});
    tmp2 = tmp2 + sum(model.Eqtkneg{k}.*(model.Eqmk2{k}+model.Eqvneg2+2*model.Eqmk{k}*model.Eqvneg));
end
model.tildeaneg = model.a00 + 0.5*tmp1;
model.tildebneg = model.b00 + 0.5*tmp2;

model.Eqtauneg = model.tildeaneg / model.tildebneg;
model.Vqtauneg = model.tildeaneg / model.tildebneg^2;
model.Eqlntauneg = psi(model.tildeaneg) - log(model.tildebneg);
%model.Hqtauneg = gammaln(model.tildeaneg) - (model.tildeaneg-1)*psi(model.tildeaneg) - log(model.tildebneg) + model.tildeaneg;

% 5. q(tau_{+1})
tmp1 = 0;
tmp2 = 0;
for k=1:K
    tmp1 = tmp1 + sum(model.Eqtkpos{k});
    tmp2 = tmp2 + sum(model.Eqtkpos{k}.*(model.Eqmk2{k}+model.Eqvpos2-2*model.Eqmk{k}*model.Eqvpos));
end
model.tildeapos = model.a00 + 0.5*tmp1;
model.tildebpos = model.b00 + 0.5*tmp2;

model.Eqtaupos = model.tildeapos / model.tildebpos;
model.Vqtaupos = model.tildeapos / model.tildebpos^2;
model.Eqlntaupos = psi(model.tildeapos) - log(model.tildebpos);
%model.Hqtaupos = gammaln(model.tildeapos) - (model.tildeapos-1)*psi(model.tildeapos) - log(model.tildebpos) + model.tildeapos;

% 6. q(tau_0)
tmp1 = 0;
tmp2 = 0;
for k=1:K
    tmp1 = tmp1 + sum(model.Eqtk0{k});
    tmp2 = tmp2 + sum(model.Eqtk0{k}.*model.Eqmk2{k});
end
model.tildea0 = model.a00 + 0.5*tmp1;
model.tildeb0 = model.b00 + 0.5*tmp2;

model.Eqtau0 = model.tildea0 / model.tildeb0;
model.Vqtau0 = model.tildea0 / model.tildeb0^2;
model.Eqlntau0 = psi(model.tildea0) - log(model.tildeb0);
%model.Hqtau0 = gammaln(model.tildea0) - (model.tildea0-1)*psi(model.tildea0) - log(model.tildeb0) + model.tildea0;

% 7. q(v_{-1})
tmp1 = 0;
tmp2 = 0;
for k=1:K
    tmp1 = tmp1 + sum(model.Eqtkneg{k});
    tmp2 = tmp2 + sum(model.Eqtkneg{k}.*model.Eqmk{k});
end

model.tildelambdapneg = model.Eqtauneg*tmp1 + model.Eqlambdaneg;
model.tildemupneg = (-model.Eqtauneg*tmp2 + model.Eqlambdaneg*model.Eqmneg)/model.tildelambdapneg;
model.tildelambdanneg = model.Eqlambdaneg;
model.tildemunneg = model.Eqmneg;

tmp3p = 0;
tmp3n = 0;
for k=1:K
    tmp3p = tmp3p + sum(model.Eqtkneg{k}.*logGauss(model.Eqmk{k}',-1,1/model.Eqtauneg)');
    tmp3n = tmp3n + sum(model.Eqtkneg{k}.*logGauss(model.Eqmk{k}',0,1/model.Eqtauneg)');
end
model.tildewpneg = tmp3p + logGauss(1,model.Eqmneg,1/model.Eqlambdaneg) - logGauss(1,model.tildemupneg,1/model.tildelambdapneg); %log(wp)
model.tildewnneg = tmp3n + logGauss(-1,model.Eqmneg,1/model.Eqlambdaneg) - logGauss(-1,model.tildemunneg,1/model.tildelambdanneg); %log(wn)
model.tildeZnegwn = 0.5*erfc(model.tildemunneg*sqrt(model.tildelambdanneg/2)) + exp(model.tildewpneg-model.tildewnneg)/2*erfc(-model.tildemupneg*sqrt(model.tildelambdapneg/2)); %Z/wn
model.tildeZnegwp = exp(model.tildewnneg-model.tildewpneg)/2*erfc(model.tildemunneg*sqrt(model.tildelambdanneg/2)) + 0.5*erfc(-model.tildemupneg*sqrt(model.tildelambdapneg/2)); %Z/wp

model.tildeMpneg0 = 1/(2*model.tildeZnegwp)*erfc(-model.tildemupneg*sqrt(model.tildelambdapneg/2));
model.tildeMpneg1 = 1/(2*model.tildeZnegwp)*(erfc(-model.tildemupneg*sqrt(model.tildelambdapneg/2))*model.tildemupneg + sqrt(2/pi/model.tildelambdapneg)/exp(model.tildelambdapneg*model.tildemupneg^2/2));
model.tildeMpneg2 = 1/(2*model.tildeZnegwp)*(erfc(-model.tildemupneg*sqrt(model.tildelambdapneg/2))*(model.tildemupneg^2+1/model.tildelambdapneg) + sqrt(2/pi/model.tildelambdapneg)*model.tildemupneg/exp(model.tildelambdapneg*model.tildemupneg^2/2));
model.tildeMnneg0 = 1/(2*model.tildeZnegwn)*erfc(model.tildemunneg*sqrt(model.tildelambdanneg/2));
model.tildeMnneg1 = 1/(2*model.tildeZnegwn)*(erfc(model.tildemunneg*sqrt(model.tildelambdanneg/2))*model.tildemunneg - sqrt(2/pi/model.tildelambdanneg)/exp(model.tildelambdanneg*model.tildemunneg^2/2));
model.tildeMnneg2 = 1/(2*model.tildeZnegwn)*(erfc(model.tildemunneg*sqrt(model.tildelambdanneg/2))*(model.tildemunneg^2+1/model.tildelambdanneg) - sqrt(2/pi/model.tildelambdanneg)*model.tildemunneg/exp(model.tildelambdanneg*model.tildemunneg^2/2));

model.Eqrneg = model.tildeMpneg1+model.tildeMnneg1;
model.Eqrneg2 = model.tildeMpneg2+model.tildeMnneg2;
model.Vqrneg = model.Eqrneg2 - model.Eqrneg^2;
model.Eqvneg = model.tildeMpneg1;
model.Eqvneg2 = model.tildeMpneg2;
model.Vqvneg = model.Eqvneg2 - model.Eqvneg^2;
%model.Hqrnegp = -((log(1/(model.tildeZnegwp*sqrt(2*pi/model.tildelambdapneg)))-0.5*model.tildemupneg^2*model.tildelambdapneg)*model.tildeMpneg0 + model.tildemupneg*model.tildelambdapneg*model.tildeMpneg1 - 0.5*model.tildelambdapneg*model.tildeMpneg2);
%model.Hqrnegn = -((log(1/(model.tildeZnegwn*sqrt(2*pi/model.tildelambdanneg)))-0.5*model.tildemunneg^2*model.tildelambdanneg)*model.tildeMnneg0 + model.tildemunneg*model.tildelambdanneg*model.tildeMnneg1 - 0.5*model.tildelambdanneg*model.tildeMnneg2);
%model.Hqrneg = model.Hqrnegp + model.Hqrnegn;

% 8. q(v_{+1})
tmp1 = 0;
tmp2 = 0;
for k=1:K
    tmp1 = tmp1 + sum(model.Eqtkpos{k});
    tmp2 = tmp2 + sum(model.Eqtkpos{k}.*model.Eqmk{k});
end

model.tildelambdappos = model.Eqtaupos*tmp1 + model.Eqlambdapos;
model.tildemuppos = (model.Eqtaupos*tmp2 + model.Eqlambdapos*model.Eqmpos)/model.tildelambdappos;
model.tildelambdanpos = model.Eqlambdapos;
model.tildemunpos = model.Eqmpos;

tmp3p = 0;
tmp3n = 0;
for k=1:K
    tmp3p = tmp3p + sum(model.Eqtkpos{k}.*logGauss(model.Eqmk{k}',1,1/model.Eqtaupos)');
    tmp3n = tmp3n + sum(model.Eqtkpos{k}.*logGauss(model.Eqmk{k}',0,1/model.Eqtaupos)');
end
model.tildewppos = tmp3p + logGauss(1,model.Eqmpos,1/model.Eqlambdapos) - logGauss(1,model.tildemuppos,1/model.tildelambdappos);
model.tildewnpos = tmp3n + logGauss(-1,model.Eqmpos,1/model.Eqlambdapos) - logGauss(-1,model.tildemunpos,1/model.tildelambdanpos);
model.tildeZposwn = 0.5*erfc(model.tildemunpos*sqrt(model.tildelambdanpos/2)) + exp(model.tildewppos-model.tildewnpos)/2*erfc(-model.tildemuppos*sqrt(model.tildelambdappos/2));
model.tildeZposwp = exp(model.tildewnpos-model.tildewppos)/2*erfc(model.tildemunpos*sqrt(model.tildelambdanpos/2)) + 0.5*erfc(-model.tildemuppos*sqrt(model.tildelambdappos/2));

model.tildeMppos0 = 1/(2*model.tildeZposwp)*erfc(-model.tildemuppos*sqrt(model.tildelambdappos/2));
model.tildeMppos1 = 1/(2*model.tildeZposwp)*(erfc(-model.tildemuppos*sqrt(model.tildelambdappos/2))*model.tildemuppos + sqrt(2/pi/model.tildelambdappos)/exp(model.tildelambdappos*model.tildemuppos^2/2));
model.tildeMppos2 = 1/(2*model.tildeZposwp)*(erfc(-model.tildemuppos*sqrt(model.tildelambdappos/2))*(model.tildemuppos^2+1/model.tildelambdappos) + sqrt(2/pi/model.tildelambdappos)*model.tildemuppos/exp(model.tildelambdappos*model.tildemuppos^2/2));
model.tildeMnpos0 = 1/(2*model.tildeZposwn)*erfc(model.tildemunpos*sqrt(model.tildelambdanpos/2));
model.tildeMnpos1 = 1/(2*model.tildeZposwn)*(erfc(model.tildemunpos*sqrt(model.tildelambdanpos/2))*model.tildemunpos - sqrt(2/pi/model.tildelambdanpos)/exp(model.tildelambdanpos*model.tildemunpos^2/2));
model.tildeMnpos2 = 1/(2*model.tildeZposwn)*(erfc(model.tildemunpos*sqrt(model.tildelambdanpos/2))*(model.tildemunpos^2+1/model.tildelambdanpos) - sqrt(2/pi/model.tildelambdanpos)*model.tildemunpos/exp(model.tildelambdanpos*model.tildemunpos^2/2));

model.Eqrpos = model.tildeMppos1+model.tildeMnpos1;
model.Eqrpos2 = model.tildeMppos2+model.tildeMnpos2;
model.Vqrpos = model.Eqrpos2 - model.Eqrpos^2;
model.Eqvpos = model.tildeMppos1;
model.Eqvpos2 = model.tildeMppos2;
model.Vqvpos = model.Eqvpos2 - model.Eqvpos^2;
%model.Hqrposp = -((log(1/(model.tildeZposwp*sqrt(2*pi/model.tildelambdappos)))-0.5*model.tildemuppos^2*model.tildelambdappos)*model.tildeMppos0 + model.tildemuppos*model.tildelambdappos*model.tildeMppos1 - 0.5*model.tildelambdappos*model.tildeMppos2);
%model.Hqrposn = -((log(1/(model.tildeZposwn*sqrt(2*pi/model.tildelambdanpos)))-0.5*model.tildemunpos^2*model.tildelambdanpos)*model.tildeMnpos0 + model.tildemunpos*model.tildelambdanpos*model.tildeMnpos1 - 0.5*model.tildelambdanpos*model.tildeMnpos2);
%model.Hqrpos = model.Hqrposp + model.Hqrposn;

% 9. q(m_{-1},lambda_{-1})
model.tildelambdamneg = 1+model.beta0;
model.tildemumneg = (model.Eqrneg+model.mu0*model.beta0)/model.tildelambdamneg;
model.tildealambdaneg = model.a01+0.5;
model.tildeblambdaneg = model.b01 + 0.5*model.beta0*model.mu0^2 + 0.5*model.Eqrneg2 - 0.5*(model.Eqrneg+model.mu0*model.beta0)^2/(model.beta0+1);

model.Eqmneg = model.tildemumneg;
model.Eqlambdaneg = model.tildealambdaneg/model.tildeblambdaneg;
model.Vqlambdaneg = model.tildealambdaneg/model.tildeblambdaneg^2;
model.Eqlnlambdaneg = psi(model.tildealambdaneg) - log(model.tildeblambdaneg);
model.Hqlambdaneg = gammaln(model.tildealambdaneg) - (model.tildealambdaneg-1)*psi(model.tildealambdaneg) - log(model.tildeblambdaneg) + model.tildealambdaneg;
%model.Hqmlambdaneg = 0.5*(1+log(2*pi)-log(model.tildelambdamneg)) - 0.5*model.Eqlnlambdaneg + model.Hqlambdaneg;

% 10. q(m_{+1}) & q(lambda_{+1})
model.tildelambdampos = 1+model.beta0;
model.tildemumpos = (model.Eqrpos+model.mu0*model.beta0)/model.tildelambdampos;
model.tildealambdapos = model.a01+0.5;
model.tildeblambdapos = model.b01 + 0.5*model.beta0*model.mu0^2 + 0.5*model.Eqrpos2 - 0.5*(model.Eqrpos+model.mu0*model.beta0)^2/(model.beta0+1);

model.Eqmpos = model.tildemumpos;
model.Eqlambdapos = model.tildealambdapos/model.tildeblambdapos;
model.Vqlambdapos = model.tildealambdapos/model.tildeblambdapos^2;
model.Eqlnlambdapos = psi(model.tildealambdapos) - log(model.tildeblambdapos);
model.Hqlambdapos = gammaln(model.tildealambdapos) - (model.tildealambdapos-1)*psi(model.tildealambdapos) - log(model.tildeblambdapos) + model.tildealambdapos;
%model.Hqmlambdapos = 0.5*(1+log(2*pi)-log(model.tildelambdampos)) - 0.5*model.Eqlnlambdapos + model.Hqlambdapos;

% 11. q(t_{j,k})
for k=1:K
    model.tilderhokneg{k} = 0.5*model.Eqlntauneg - 0.5*model.Eqtauneg*(model.Eqmk2{k} + model.Eqvneg2 + 2*model.Eqvneg*model.Eqmk{k}) - log(2);
    model.tilderhokpos{k} = 0.5*model.Eqlntaupos - 0.5*model.Eqtaupos*(model.Eqmk2{k} + model.Eqvpos2 - 2*model.Eqvpos*model.Eqmk{k}) - log(2);
    model.tilderhok0{k} = 0.5*model.Eqlntau0 - 0.5*model.Eqtau0*model.Eqmk2{k};
    logdiff = -S{k}*model.Eqw; %J_k x 1

    model.Eqtkneg{k} = (exp(model.tilderhokneg{k}) ./ (exp(model.tilderhokneg{k}) + exp(model.tilderhokpos{k}) + exp(model.tilderhok0{k}+logdiff))) .* (sum(S{k}(:,1:1),2)~=0);
    model.Eqtkpos{k} = (exp(model.tilderhokpos{k}) ./ (exp(model.tilderhokneg{k}) + exp(model.tilderhokpos{k}) + exp(model.tilderhok0{k}+logdiff))) .* (sum(S{k}(:,1:1),2)~=0);
    %model.Eqtkneg{k} = (exp(model.tilderhokneg{k}) ./ (exp(model.tilderhokneg{k}) + exp(model.tilderhokpos{k}) + exp(model.tilderhok0{k}+logdiff)));
    %model.Eqtkpos{k} = (exp(model.tilderhokpos{k}) ./ (exp(model.tilderhokneg{k}) + exp(model.tilderhokpos{k}) + exp(model.tilderhok0{k}+logdiff)));
    model.Eqtk0{k} = 1-model.Eqtkneg{k}-model.Eqtkpos{k};
    %model.Hqtk{k} = -(model.Eqtkneg{k}.*log(model.Eqtkneg{k}+1e-6) + model.Eqtkpos{k}.*log(model.Eqtkpos{k}+1e-6) + model.Eqtk0{k}.*log(model.Eqtk0{k}+1e-6));
end

% 12. q(w)
tmp1 = 0;
tmp2 = 0;
for k=1:K
    tmp1 = tmp1 + bsxfun(@times, fchi(model.xi{k}), S{k})'*S{k}; %D x D
    tmp2 = tmp2 + sum(bsxfun(@times, model.Eqtkneg{k}+model.Eqtkpos{k}-0.5, S{k}), 1)'; %D x 1
end

model.tildeLambdaw = model.EqLambda + 2*tmp1;
U = chol(model.tildeLambdaw);
model.tildemuw = U \ (U' \ tmp2);

model.Eqw = model.tildemuw;
model.Vqw = U \ (U' \ eye(D));
model.Eqw2 = model.Vqw + model.Eqw*model.Eqw';
%model.Hqw = -sum(log(diag(U))) + 0.5*D*(1+log(2*pi));

% 13. q(Lambda)
model.tildenuLambda = model.nu0 + 1;
model.tildeWLambdainv = model.W0inv + model.Eqw2;

U = chol(model.tildeWLambdainv);
model.tildeWLambda = U \ (U' \ eye(D));

model.EqLambda = model.tildenuLambda * model.tildeWLambda;
model.EqlnLambda = -2*sum(log(diag(U))) + D*log(2) + sum(psi((model.tildenuLambda+1-(1:D))/2));
%logB = model.tildenuLambda*sum(log(diag(U))) - model.tildenuLambda*D/2*log(2) + D*(D-1)/4*log(pi) + sum(gammaln((model.tildenuLambda+1-(1:D))/2));
%model.HqLambda = -logB - (model.tildenuLambda-D-1)/2*model.EqlnLambda + model.tildenuLambda*D/2;

end

function model = maximize(S, K, model)
%% M-step
for k=1:K
    model.xi{k} = sqrt(diag(S{k}*model.Eqw2*S{k}')); %J_k x 1
end

end

function llh = eval(S, Z, Ld, Ik, K, N, model)
%% Evaluation for algorithm termination
% compute the expected complete-data log likelihood E[logP(x,z)]
D = size(S{1},2);

%1. p(z|u)
tmp1 = 0;
for k=1:K
    tmp1 = tmp1 + sqrt(N)*Z{k}'*model.Equk{k} - 0.5*N*trace(model.Equk2{k}*Ld{k}) - 0.5*model.lnLd(k);
end

%2. p(u|m,lambda)
tmp2 = 0;
for k=1:K
    zs = diag(model.Equk2{k});
    for j=1:length(Ik{k})
        tmp2 = tmp2 + sum(model.Eqlnlambdak{k}(j) - 0.5*model.Eqlambdak{k}(j) * ( zs(sum(Ik{k}(1:j))-Ik{k}(j)+1 : sum(Ik{k}(1:j))) ...
            +model.Eqmk2{k}(j)-2*model.Eqmk{k}(j)*model.Equk{k}(sum(Ik{k}(1:j))-Ik{k}(j)+1 : sum(Ik{k}(1:j))) ));
    end
end

%3. p(lambda)
tmp3 = 0;
for k=1:K
    tmp3 = tmp3 + sum((model.a00-1)*model.Eqlnlambdak{k}-model.b00*model.Eqlambdak{k});
end

%4. p(m|t,r,tau)
tmp4=0;
for k=1:K
    tmp4 = tmp4 + sum( ...
    model.Eqtkneg{k}.*(-0.5*model.Eqtauneg*(model.Eqmk2{k}+model.Eqvneg2+2*model.Eqvneg*model.Eqmk{k})+0.5*model.Eqlntauneg-0.5*log(2*pi)) + ...
    model.Eqtk0{k}.*(-0.5*model.Eqtau0*model.Eqmk2{k}+0.5*model.Eqlntau0-0.5*log(2*pi)) + ...
    model.Eqtkpos{k}.*(-0.5*model.Eqtaupos*(model.Eqmk2{k}+model.Eqvpos2-2*model.Eqvpos*model.Eqmk{k})+0.5*model.Eqlntaupos-0.5*log(2*pi)) ...
    );
end

%5. p(tau)
tmp5 = (model.a00-1)*model.Eqlntauneg - model.b00*model.Eqtauneg + ...
    (model.a00-1)*model.Eqlntau0 - model.b00*model.Eqtau0 + ...
    (model.a00-1)*model.Eqlntaupos - model.b00*model.Eqtaupos;

%6. p(r|m,lambda)
tmp6 = 0.5*model.Eqlnlambdaneg - 0.5*model.Eqlambdaneg*(model.Eqrneg2-2*model.Eqrneg*model.Eqmneg) - 0.5*model.Eqlambdaneg*model.Eqmneg^2 - 0.5/model.tildelambdamneg + ...
    0.5*model.Eqlnlambdapos - 0.5*model.Eqlambdapos*(model.Eqrpos2-2*model.Eqrpos*model.Eqmpos) - 0.5*model.Eqlambdapos*model.Eqmpos^2 - 0.5/model.tildelambdampos;

%7. p(m,lambda)
tmp7 = -0.5*model.beta0/model.tildelambdamneg - 0.5*model.beta0*model.Eqlambdaneg*(model.Eqmneg^2 + model.mu0^2 - 2*model.mu0*model.Eqmneg) + ...
    0.5*model.Eqlnlambdaneg - model.b01*model.Eqlambdaneg + (model.a01-1)*model.Eqlnlambdaneg ...
    -0.5*model.beta0/model.tildelambdampos - 0.5*model.beta0*model.Eqlambdapos*(model.Eqmpos^2 + model.mu0^2 - 2*model.mu0*model.Eqmpos) + ...
    0.5*model.Eqlnlambdapos - model.b01*model.Eqlambdapos + (model.a01-1)*model.Eqlnlambdapos;
    
%8. p(t|w)
tmp8 = 0;
for k=1:K
    tmp8 = tmp8 + sum(...
        (model.Eqtkneg{k}+model.Eqtkpos{k})*log(0.5) + ...
        S{k}*model.Eqw.*(model.Eqtkneg{k}+model.Eqtkpos{k}) + ...
        log(sigmoid(model.xi{k})) + ...
        -0.5*(S{k}*model.Eqw+model.xi{k}) - fchi(model.xi{k}).*(diag(S{k}*model.Eqw2*S{k}')-model.xi{k}.^2) ...
        );
end

%9. p(w|Lambda)
tmp9 = -0.5*trace(model.Eqw2*model.EqLambda) + 0.5*model.EqlnLambda;

%10. p(Lambda)
tmp10 = -0.5*trace(model.W0inv*model.EqLambda) + 0.5*(model.nu0-D-1)*model.EqlnLambda;

llh = tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8+tmp9+tmp10;

end




