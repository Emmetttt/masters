lowZred = [];
lowZdist = [];
lowZRA = [];
lowZdec = [];
lowZdisterror = [];
lowZmpc = [];
lowZmpcerror = [];
lowZH = [];


for i=1:length(redshift)
    if redshift(i) < 0.1
        if redshift(i) > 0.03
            lowZred = [lowZred redshift(i)];
            lowZdist = [lowZdist distmod(i)];
            lowZRA = [lowZRA ra(i)];
            lowZdec = [lowZdec dec(i)];
            lowZdisterror = [lowZdisterror distmoderror(i)];
            lowZmpc = [lowZmpc distmpc(i)];
            lowZmpcerror = [lowZmpcerror distmpcerror(i)];
            lowZH = [lowZH H(i)];
        end
    end
end

H0 = 70;
deltaH = [];

for j=1:length(lowZred)
    deltaH = [deltaH H0-lowZH(j)];
end