function [ F ] = kibam_c_k(x)

global battery;

c(1) = x(1);
k(1) = x(2);


F(1) = c(1) - (battery.F1*(1-exp(-k(1)*battery.t)*battery.t1) - (1-exp(-k(1)*battery.t1)*battery.t) ...
                / (battery.F1*(1-exp(-k(1)*battery.t)*battery.t1) - (1-exp(-k(1)*battery.t1)*battery.t) - ...
                k(1)*battery.F1*battery.t*battery.t1 + k(1)*battery.t*battery.t1));


F(2) = c(1) - (battery.F2*(1-exp(-k(1)*battery.t)*battery.t2) - (1-exp(-k(1)*battery.t2)*battery.t) ...
                / (battery.F2*(1-exp(-k(1)*battery.t)*battery.t2) - (1-exp(-k(1)*battery.t2)*battery.t) - ...
                k(1)*battery.F2*battery.t*battery.t2 + k(1)*battery.t*battery.t2));   



end

