

function rate = miniMap(bitPat)
%input is the bit pattern that determines the modulation and coding scheme
%as defined by 802.11-2012, returns the corresponding data rate.
%function written because I didnt think ahead
if isequal(bitPat,[1 1 0 1])
    rate = 6;
elseif isequal(bitPat,[1 1 1 1])
    rate = 9;
elseif isequal(bitPat,[0 1 0 1])
    rate = 12;
elseif isequal(bitPat,[0 1 1 1])
    rate = 18;
elseif isequal(bitPat,[1 0 0 1])
    rate = 24;
elseif isequal(bitPat,[1 0 1 1])
    rate = 36;
elseif isequal(bitPat,[0 0 0 1])
    rate = 48;
elseif isequal(bitPat,[0 0 1 1])
    rate = 54;
end
end