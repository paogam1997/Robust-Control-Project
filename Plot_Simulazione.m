%% plot simulazione
teta1=out.teta1.Data ;
teta2=out.teta2.data;
% t=out.teta1.Time;
figure(); hold on;grid on; xlim([-1,1]); ylim([-1,1]);
p0=[0;0];
plot_origine = plot(p0(1),p0(2),'or');
N=size(teta2,1);
for k=1:N
    t=10*k/N;
    q=[teta1(k),teta2(k)];
    [p1,p2]=cinematica(p0,q(1),q(2),l1,l2);
    plot_link1=plot([p0(1),p1(1)],[p0(2),p1(2)],'b');
    plot_giunto2=plot(p1(1),p1(2),'or');
    plot_link2=plot([p1(1),p2(1)],[p1(2),p2(2)],'b');
    plot_ee=plot(p2(1),p2(2),'*k');
    timer=text(-1,1,"Timer: "+num2str(t,2));
    pause(0.01);
    delete(plot_link1);
    delete(plot_giunto2);
    delete(plot_link2);
    delete(plot_ee);
    delete(timer);
end
function [p1,p2]=cinematica(p0, q1, q2,l1,l2)
% q1=teta1;
% q2=teta2;
p1 = [l1*sin(q1); l1*cos(q1)];
p2 = p0 + p1 + [l2*sin(q2); l2*cos(q2)];
end