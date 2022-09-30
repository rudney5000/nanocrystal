function Ei = Field(T, Tm, Am, t)
if and((0 <= t), (t <= T))
   Ei = Am*t;
   else
   Ei = Am*Tm - Am*t;
end
end