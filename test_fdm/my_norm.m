function [out1,out2] = my_norm(in1,in2)
%MY_NORM Normalizes in1 and in2.

%in1 = abs(in1);
%in2 = abs(in2);

%max_fac = max(in1,[],1);
max_fac = in1(1,:);
out1 = bsxfun(@rdivide,in1,max_fac);
out2 = bsxfun(@rdivide,in2,max_fac);



end

