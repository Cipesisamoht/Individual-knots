poly = Jones([1,1,1])

function x = Kaufmann(braid)
pos = abs(braid);
dir = braid./pos;
len = length(braid);
checkarray = 1:max(pos);
poly1 = [];
for i = 1:len
    if dir(i) == 1
        poly1 = [poly1, 0, pos(i)];
    else
        poly1 = [poly1, pos(i), 0];
    end
end
loops1 = 2^len;
loopcache = [];
ucache = [];
for i = 1:loops1
    k = i;
    for j = 1:len
        if k>loops1/(2^j)
            ucache = [ucache, poly1(2*j-1)];
            k=k-loops1/(2^j);
        else
            ucache = [ucache, poly1(2*j)];
        end
    end
    loops = 0;
    ucache = ucache(ismember(ucache,setdiff(ucache,[0])));
    ucache = sort(ucache);
    loops = loops+length(setdiff(checkarray,ucache));
    if length(ucache)>1
        for j = 1:length(ucache)-1
            first = ucache(j);
            second = ucache(j+1);
            if second-first ~= 1
                loops = loops + 1;
            end
        end
    end
    loopcache = [loopcache, loops];
    ucache = [];
end
poly = [0]; %decode loopcache
first = 0;
last = 0;
for i = 1:length(loopcache)
    position = 0;
    currentI = i;
    for j = 1:len
        floatLength = length(loopcache);
        if currentI > floatLength/(2^j)
            position = position + 1;
            currentI = currentI-floatLength/(2^j);
        else
            position = position -1;
        end
    end
    positions = [position]; %finds the starting position of each elment of loopCache
    positionsNew = [];
    multiply = 1;
    numberof = 1;
    if loopcache(i) > 0
        for j = 1:loopcache(i)
            multiply = -multiply;
            for k = 1:numberof
                positionsNew = [positionsNew, positions(k)-2, positions(k)+2];
                positionsNew = unique(positionsNew);
            end
            numberof = numberof+1;
            positions = positionsNew;
            positionsNew = [];
        end
    else
        positions = [position];
    end %starts generating the polynomial
    oldFirst = first;
    first = min([first, positions]);
    if oldFirst ~= first
        difference = oldFirst-first;
        poly = [zeros(1,difference), poly];
    end
    oldLast = last;
    last = max([last,positions]);
    if oldLast ~= last
        difference = last-oldLast;
        poly = [poly,zeros(1,difference)];
    end
    for j = 1:length(positions)
        current = positions(j);
        poly(current-first+1) = poly(current-first+1)+multiply*nchoosek(length(positions)-1,j-1); %writes into the poly function
    end
end
x = [-first, -poly];
end

function x = Jones(braid)% uses the kauffman function to find the Jones polynomial. output form is [highest power of A, coeficient of highest power, coeficient of second highest power ...]
poly=Kaufmann(braid);
writhe1 = -sum(braid./abs(braid));
order = poly(1);
poly = poly(2:end);
poly = poly.*((-1)^writhe1);
order = order-3*writhe1;
x = [order,poly];
end