function [reqs, xeqs, rcs, xms] = tptm(v_ratio, config, vl_oc, il_oc, p3_oc, vl_sc, il_sc, p3_sc)

% v_ratio is line primary voltage divided by line secondary voltage

pph_oc = p3_oc / 3;             % per phase OC power
pph_sc = p3_sc / 3;             % per phase SC power

if config == 1                          % Y-Y

        a = v_ratio                     % per phase a factor
        
        vph_oc = vl_oc / sqrt(3);       % open circuit phase voltage
        iph_oc = il_oc;                 % open circuit phase current
        
        vph_sc = vl_sc / sqrt(3);       % short circuit phase voltage
        iph_sc = il_sc;                 % short circuit phase current

else if config == 2                     % Y-D

        a = v_ratio / sqrt(3)           % per phase a factor
        
        if v_ratio < 1                  % OC test performed on the primary (Y) side
            
            vph_oc = vl_oc / sqrt(3);   % open circuit phase voltage
            iph_oc = il_oc;             % open circuit phase current
            
            % and SC test performed on secondary (D) side
            vph_sc = vl_sc;             % short circuit phase voltage
            iph_sc = il_sc / sqrt(3);   % short circuit phase current

        else if v_ratio > 1             % OC test performed on secondary (D) side
           
            vph_oc = vl_oc;             % open circuit phase voltage
            iph_oc = il_oc / sqrt(3);   % open circuit phase current
            
            % and SC test performed on primary (Y) side

            vph_sc = vl_sc / sqrt(3);   % short circuit phase voltage
            iph_sc = il_sc;             % short circuit phase current

        end

else if config == 3                     % D-Y

        a = v_ratio * sqrt(3)           % per phase a factor
        
         if v_ratio < 1                 % OC test performed on the primary (D) side
            
            vph_oc = vl_oc;             % open circuit phase voltage
            iph_oc = il_oc / sqrt(3);   % open circuit phase current
            
            % and SC test performed on secondary (Y) side

            vph_sc = vl_sc / sqrt(3);   % short circuit phase voltage
            iph_sc = il_sc;             % short circuit phase current

        else if v_ratio > 1             % OC test performed on secondary (Y) side
           
            vph_oc = vl_oc / sqrt(3);   % open circuit phase voltage
            iph_oc = il_oc;             % open circuit phase current
            
            % and SC test performed on primary (D) side

            vph_sc = vl_sc;             % short circuit phase voltage
            iph_sc = il_sc / sqrt(3);   % short circuit phase current
       
        end


         else                           % D-D, inputting 4 as config is exclusive to 1,2, & 3

        a = v_ratio                     % per phase a factor
        
        vph_oc = vl_oc;                 % open circuit phase voltage
        iph_oc = il_oc / sqrt(3);       % open circuit phase current
        
        vph_sc = vl_sc;                 % short circuit phase voltage
        iph_sc = il_sc / sqrt(3);       % short circuit phase current

end

Y_shunt = (iph_oc / vph_oc) * exp(-j*acos(pph_oc/(vph_oc*iph_oc)));
rc = 1/real(Y_shunt);
xm = -1/imag(Y_shunt);

Z_eq = (vph_sc / iph_sc) * exp(j*acos(pph_sc / (vph_sc*iph_sc)));
req = real(Z_eq);
xeq = imag(Z_eq);

if v_ratio < 1                          % OC performed on primary, SC performed on secondary, relative to 2ndary model
    rcs = rc / (a^2);
    xms = xm / (a^2);
    reqs = req;
    xeqs = xeq;

else if v_ratio < 1                     % OC performed on secondary, SC performed on primary, relative to 2ndary model
    rcs = rc;
    xms = xm;
    reqs = req / (a^2);
    xeqs = xeq / (a^2);
% config = [1, 2, 3, 4]; % Y-Y, Y-D, D-Y, D-D

