function [k, E] = getTurbSpec(L)

uvel = reshape(L(1:512*512, :), 512, 512, 1000);
vvel = reshape(L(512*512+1:2*512*512, :), 512, 512, 1000);
wvel = reshape(L(2*512*512+1:end, :), 512, 512, 1000);

[k, E] = turbspec(uvel, vvel, wvel);
end