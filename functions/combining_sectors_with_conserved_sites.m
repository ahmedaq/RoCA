function s_out = combining_sectors_with_conserved_sites(s_in,cs,nns,ls)

% Code to include all neighboring conserved sites in the sectors
% 
% Last updated: 2017-03-25
% 
% INPUT
%     s_in                  Sector
%     cs                    Fully conserved sites
%     nns                   Parameter to define neighborhood. For example:
%                           - nns = 1 --> check only the immediate left
%                           and right neighbor of a sector site. If sector
%                           site = 5, neighboring sites = 4 and 6.
%                           - nns = 2 --> check two neighboring sites on
%                           each side of a sector site. If sector site = 5,
%                           neiboring sites = 3, 4, 6, and 7.
%     ls                    Correct number of sites in the protein (before
%                           preprocessing)
%     
% OUTPUT
%     s_out                 Sector supplemented with neighboring conserved
%                           sites

%%

lsin = length(s_in);
s_out = [];
for k = 1:lsin
    ns_left = max(s_in(k)-nns,1):1:s_in(k)-1;
    ns_c_left = cs(ismember(cs,ns_left)); %neighboring sites that are conserved
    ns_c_new_left2 = [];
    for mm = 1:length(ns_c_left)
        test_cs = ns_c_left(mm);
        ns_c_new_left = [];
        left = 1;
        while left == 1
            test_cs_left = max(test_cs,1)-1;
            if ismember(test_cs_left,cs)==1
                ns_c_new_left = [ns_c_new_left test_cs_left];
                test_cs = test_cs - 1;
            else
                left = 2;
            end
        end
        ns_c_new_left2 = [ns_c_new_left2 ns_c_new_left];
    end
    s_out = [s_out ns_c_left unique(ns_c_new_left2)];
    
    ns_right = s_in(k)+1:1:min(s_in(k)+nns,ls);
    ns_c_right = cs(ismember(cs,ns_right)); %neighboring sites that are conserved
    ns_c_new_right2 = [];
    for mm = 1:length(ns_c_right)
        test_cs = ns_c_right(mm);
        ns_c_new_right = [];
        right = 1;
        while right == 1
            test_cs_right = min(test_cs,ls)+1;
            if ismember(test_cs_right,cs)==1
                ns_c_new_right = [ns_c_new_right test_cs_right];
                test_cs = test_cs + 1;
            else
                right = 2;
            end
        end
        ns_c_new_right2 = [ns_c_new_right2 ns_c_new_right];
    end
    s_out = [s_out ns_c_right unique(ns_c_new_right2)];
    
    s_out = [unique([s_out s_in(k)])];
    
end     