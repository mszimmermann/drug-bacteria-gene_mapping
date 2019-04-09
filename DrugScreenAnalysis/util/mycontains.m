function containflag = contains(x, y)
containflag = ~isempty(strfind(x,y));
%if query is empty, subject contains it
if isempty(y)
    containflag = (1==1);
end
