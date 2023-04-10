function angle = pi2pi(angle)

n = length(angle);

for i=1:n
    if (angle(i)<-2*pi) || (angle(i)>2*pi)
        angle(i) = mod(angle(i), 2*pi);
    end
    
    if (angle(i)>pi)
        angle(i) = angle(i)-2*pi;
    end
    
    if (angle(i)<-pi)
        angle(i) = angle(i)+2*pi;
    end
    
end

end