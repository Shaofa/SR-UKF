function [ y ] = srukf_h_func(x)
    y(1,:) = sqrt(x(1,:).^2+x(4,:).^2);
    y(2,:) = atan2(x(4,:),x(1,:));
    y(3,:) = (x(1,:).*x(2,:)+x(4,:).*x(5,:))./y(1,:);
end