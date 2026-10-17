clear; clc
%Within-host model of bacterial meningitis in the cerebrospinal fluid (CSF):
%B(t) = bacterial concentration (CFU/mL)
%I(t) = immune response, measured as the CSF interleukin-6 (IL-6) concentration (pg/mL)
%
%dB/dt = r*B*(1 - B/K) - k*B*I               (1)  logistic growth minus immune killing
%dI/dt = s*B/(h + B) - d*(I - I0)            (2)  saturating IL-6 production minus
%                                                      decay back to the baseline level I0
%
% What this script does
%     1. Sets the parameters (units and sources are given next to each value).
%     2. Calibrates the immune killing rate k and the maximum IL-6 production rate s.
%     3. Computes the threshold quantity R0 and the two equilibria.
%     4. Checks that both equilibria satisfy the model equations.
%     5. Computes the Jacobian and its eigenvalues at each equilibrium (local stability).
%     6. Tests how the stability of the persistent equilibrium depends on the two
%        assumed parameters h and d.

%%\ Parameters
r  = 15;            % bacterial growth rate (per day)
K  = 1e9;           % bacterial carrying capacity (CFU/mL)
d  = 0.5;           % per day (placeholder)
I0 = 205.83;        % baseline CSF IL-6 (pg/mL)
Is = 3538.69;       % CSF IL-6 at the persistent equilibrium (pg/mL)
Bs = 1e6;           % bacterial load at the persistent equilibrium (CFU/mL) 
h  = 1e6;           % half-saturation constant (CFU/mL)
d = 0.5;            % IL-6 decay rate (per day)

%%\ Calibrated parameters
k = r*(1 - Bs/K)/Is               % immune killing rate (mL/(pg*day)) %[output:43dc6465]
s = d*(Is - I0)*(h + Bs)/Bs        % maximum IL-6 production rate (pg/(mL day)) %[output:421949aa]
%%\ Threshold
R0 = r/(k*I0) %[output:54a34372]
%%\ Disease-free equilibrium 
B_dfe = 0;
I_dfe = I0;

%%\ chronic equilibrium
a  = r/K;
b  = -(r*(1 - h/K) - k*I0 - k*s/d);
c  = -(r - k*I0)*h;
sg = 2*(b >= 0) - 1;
q  = -(b + sg*sqrt(b^2 - 4*a*c))/2;      
B_ee = max(q/a, c/q);         % the positive root (c < 0 when r > k*I0)
I_ee = I0 + (s/d)*B_ee/(h + B_ee); % chronic equilibrium exists iff R0 > 1

fprintf('DFE: (B, I) = (%.4g, %.4f)\n', B_dfe, I_dfe); %[output:557e1572]
if r > k*I0 %[output:group:36cffe34]
    fprintf('EE:  (B, I) = (%.6g, %.4f)\n', B_ee, I_ee); %[output:7889fd75]
else
    fprintf('No persistent equilibrium (r <= k*I0)\n');
end %[output:group:36cffe34]
%%\ Check
f = @(B,I) r*B*(1 - B/K) - k*B*I;
g = @(B,I) s*B/(h + B) - d*(I - I0);
check_dfe = [f(B_dfe,I_dfe), g(B_dfe,I_dfe)] %[output:85b8bab9]
check_ee  = [f(B_ee,I_ee),   g(B_ee,I_ee)] %[output:6b73ee6d]
%%\ Jacobian
% Near an equilibrium, small disturbances grow or shrink according to the eigenvalues of
% the Jacobian matrix. If all eigenvalues have negative real part, the equilibrium is
% locally stable. Complex eigenvalues mean they return to equilibrium oscillates (spiral).
J = @(B,I) [r*(1 - 2*B/K) - k*I, -k*B;
    s*h/(h + B)^2,       -d];

J_dfe = J(B_dfe, I_dfe) %[output:68d091ff]
J_ee  = J(B_ee, I_ee) %[output:17cfb65b]
%%\ Eigenvalues
eig_dfe = eig(J_dfe) %[output:2edd2bf2]
eig_ee  = eig(J_ee) %[output:51a338f7]
%%\ Analytic check
eig_dfe_formula = [r - k*I0; -d] %[output:6ca092e7]
trace_ee = -r*B_ee/K - d %[output:068bce09]
det_ee   = (r*B_ee/K)*d + k*B_ee*s*h/(h + B_ee)^2 %[output:3e481f54]
[trace(J_ee), det(J_ee)] %[output:4b172dba]
%%\ Stability
if all(real(eig_dfe) < 0), disp('DFE: stable'), else, disp('DFE: unstable'), end %[output:6235f0bc]
if all(real(eig_ee) < 0),  disp('EE: stable'),  else, disp('EE: unstable'),  end %[output:4ed975ad]
if any(imag(eig_ee) ~= 0), disp('EE: spiral'), else, disp('EE: node'), end %[output:052373b4]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":6.6}
%---
%[output:43dc6465]
%   data: {"dataType":"textualVariable","outputData":{"name":"k","value":"0.0042"}}
%---
%[output:421949aa]
%   data: {"dataType":"textualVariable","outputData":{"name":"s","value":"3.3329e+03"}}
%---
%[output:54a34372]
%   data: {"dataType":"textualVariable","outputData":{"name":"R0","value":"17.2095"}}
%---
%[output:557e1572]
%   data: {"dataType":"text","outputData":{"text":"DFE: (B, I) = (0, 205.8300)\n","truncated":false}}
%---
%[output:7889fd75]
%   data: {"dataType":"text","outputData":{"text":"EE:  (B, I) = (1e+06, 3538.6900)\n","truncated":false}}
%---
%[output:85b8bab9]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"check_dfe","rows":1,"type":"double","value":[["0","0"]]}}
%---
%[output:6b73ee6d]
%   data: {"dataType":"matrix","outputData":{"columns":2,"exponent":"-8","name":"check_ee","rows":1,"type":"double","value":[["0.1863","0"]]}}
%---
%[output:68d091ff]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"J_dfe","rows":2,"type":"double","value":[["14.1284","0"],["0.0033","-0.5000"]]}}
%---
%[output:17cfb65b]
%   data: {"dataType":"matrix","outputData":{"columns":2,"exponent":"3","name":"J_ee","rows":2,"type":"double","value":[["-0.0000","-4.2346"],["0.0000","-0.0005"]]}}
%---
%[output:2edd2bf2]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"eig_dfe","rows":2,"type":"double","value":[["-0.5000"],["14.1284"]]}}
%---
%[output:51a338f7]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"eig_ee","rows":2,"type":"complex","value":[["-0.2575 + 1.8627i"],["-0.2575 - 1.8627i"]]}}
%---
%[output:6ca092e7]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"eig_dfe_formula","rows":2,"type":"double","value":[["14.1284"],["-0.5000"]]}}
%---
%[output:068bce09]
%   data: {"dataType":"textualVariable","outputData":{"name":"trace_ee","value":"-0.5150"}}
%---
%[output:3e481f54]
%   data: {"dataType":"textualVariable","outputData":{"name":"det_ee","value":"3.5358"}}
%---
%[output:4b172dba]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["-0.5150","3.5358"]]}}
%---
%[output:6235f0bc]
%   data: {"dataType":"text","outputData":{"text":"DFE: unstable\n","truncated":false}}
%---
%[output:4ed975ad]
%   data: {"dataType":"text","outputData":{"text":"EE: stable\n","truncated":false}}
%---
%[output:052373b4]
%   data: {"dataType":"text","outputData":{"text":"EE: spiral\n","truncated":false}}
%---
