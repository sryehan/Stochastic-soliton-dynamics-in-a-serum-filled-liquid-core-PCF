function s = ifelse_str(c1, s1, c2, s2, s3)
%IFELSE_STR  Small helper: returns s1 if c1, else s2 if c2, else s3.
if c1, s = s1; elseif c2, s = s2; else, s = s3; end
end
