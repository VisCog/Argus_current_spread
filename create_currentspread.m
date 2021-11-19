function [I, p, Vq] = create_currentspread(p)

if ~isfield(p, 'a1')
    a1 = p.a; a2 = p.a;
else
    a1 = p.a1; a2 = p.a2;
end

CM = 1./(1+((p.k*p.R1).^p.gamma));
I1 = a1.*CM;
CM = 1./(1+((p.k*p.R2).^p.gamma));
I2 = a2.*CM;

I = I1+I2;

if isfield(p, 'loc')
    Vq = interpn(unique(p.Y), unique(p.X), unique(p.Z), I, p.loc(2), p.loc(1), p.loc(3));
else
    Vq = NaN;
end
end