%将"用户(user)”改成了“车辆(vehicle)”
%为车辆添加了速度属性speed,其取值范围为[60,100] km/h。速度会影响车辆的本地计算时延和传输时延。
%用Matlab实现了修改后的模型和算法PTAS-TAHRC。主要修改包括:
%1初始化车辆数量,速度,本地计算能力等参数
%2根据速度计算本地计算时延和传输时延
%3调用PTAS-TAHRC进行资源分配和卸载决策
%4输出结果,包括车辆分配结果、能耗等

% 参数设置
num_vehicles = 50;  
speed_range = [60,100];
local_capability_range = [2 4]*1e9;  
server_capability = 20*1e9;
memory_size = 100*1e6;
required_cycles_range = [1e8 1e9];
data_size_range = [1e6 1e7];
deadline_range = [100 200];

% 车辆参数初始化
speed = unidrnd(speed_range(2)-speed_range(1)+1,num_vehicles,1) + speed_range(1) - 1; 
local_capability = unidrnd(local_capability_range(2)-local_capability_range(1)+1,num_vehicles,1) + local_capability_range(1) - 1;
required_cycles = unidrnd(required_cycles_range(2)-required_cycles_range(1)+1,num_vehicles,1) + required_cycles_range(1) - 1;
data_size = unidrnd(data_size_range(2)-data_size_range(1)+1,num_vehicles,1) + data_size_range(1) - 1;
deadline = unidrnd(deadline_range(2)-deadline_range(1)+1,num_vehicles,1)./1000 + deadline_range(1)/1000 - .01;

% 计算时延
B = 20*10^6; % 带宽 20MHz
rho = 0.5; % 传输功率 0.5W
sigma = 2*10^-13; % 噪声功率
location = rand(num_vehicles,2); 
distance = sqrt(sum(location.^2));
h = channel_gain(distance); % 根据距离计算信道增益

r = B*log2(1 + rho.*h./sigma.^2); % 计算上行传输速率
transmission_delay = data_size./r;
local_delay = required_cycles./local_capability;

% 能耗参数    
transmission_energy_coefficient = 1e-28; 
computation_energy_coefficient = 1e-26;

% 调用任务分配算法
[allocation, total_energy] = PTASTAHRC(num_vehicles, speed,...
                                   local_capability, transmission_delay,...
                                   required_cycles, data_size,... 
                                   deadline, server_capability,...
                                   memory_size,...
                                   transmission_energy_coefficient,...
                                   computation_energy_coefficient,local_delay);
                           
% 输出结果
disp('Allocation results:');
disp(allocation);
disp('Total energy consumption:');  
disp(total_energy);
% 生成结果表
table = table(speed,allocation,required_cycles);

% 绘制任务量分布
figure;
histogram(table.required_cycles)

% 绘制远程执行比例  
pie(sum(table.allocation)/50) 

% 绘制速度与执行位置的关系
figure;
scatter(table.speed, table.allocation);

% 汇总能耗
energy = [sum(E_l.*alloc) sum(E_r.*(1-alloc)) sum(E_l+E_r)];
% 任务分配算法
function [allocation, total_energy] = PTASTAHRC(num_vehicles, speed, local_capability,...
    transmission_delay, required_cycles, data_size, deadline, server_capability,...
    memory_size, transmission_energy_coefficient, computation_energy_coefficient,...
    local_delay)

    % 初始化
    N = 1:num_vehicles; 
    K = 20; % 服务器可支持的最大上传链接数
    c = [server_capability, memory_size]; % 服务器计算和存储资源
    
    E_l = computation_energy_coefficient*local_capability.^2.*required_cycles; % 本地计算能耗
    E_r = transmission_energy_coefficient.*data_size./speed; % 传输能耗
    
    % 资源分配
    d_i1 = zeros(num_vehicles, 1);
    for i = 1:num_vehicles
      d_i1(i)= required_cycles(i)/(deadline(i) - transmission_delay(i)); 
    end
    allocation = zeros(num_vehicles, 1);
    [E_l, index] = sort(E_l,'descend'); 
    for i = 1:50
      
      if all(local_delay(i) > deadline(i)) && K > 0 && all(c >= [d_i1(i),data_size(i)])
        allocation(i) = 1;  
        K = K - 1;
        c = c - [d_i1(i),data_size(i)];
      end
    end
    
    % 卸载决策
    epsilon = 0.5; 
    q = floor(2/epsilon);
    N_r = [];
    V = 0;
    for Q = nchoosek(1:num_vehicles,q)
        if length(Q) == 1
          if d_i1(Q) <= server_capability && data_size(Q) <= memory_size
            S = setdiff(1:num_vehicles,Q); 
            V_temp = sum(E_s(Q)) + GREEDYALLOC(S, E_s(S), c-sum([d_i1(Q) data_size(Q)]), K-length(Q));
            if V_temp > V
              V = V_temp;
              N_r = [Q; S(GREEDYALLOC(S,E_s(S)))]; 
            end            
          end
        else
          if sum(d_i1(Q)) <= server_capability && sum(data_size(Q)) <= memory_size  
            S = setdiff(1:num_vehicles,Q); 
            V_temp = sum(E_s(Q)) + GREEDYALLOC(S, E_s(S), c-sum([d_i1(Q) data_size(Q)]), K-length(Q));
            if V_temp > V
              V = V_temp;
              N_r = [Q; S(GREEDYALLOC(S,E_s(S)))]; 
            end
          end
        end
    end
    
    % 返回结果
    nrs= sort(N_r);
    allocation(nrs) = 1;  
    total_energy = sum(E_l.*allocation) + sum(E_r.*(1-allocation));

end

% 贪心分配子函数
function [E_s, index] = GREEDYALLOC(S, E_s, c, K)
  
  % LP松弛解
  %这部分利用MATLAB的线性规划函数linprog来获得线性规划松弛后的解。
    %linprog的输入参数为空,表示自由线性规划,没有限制条件
    %仅有一个目标函数系数c,对应GREEDY-ALLOC中传入的剩余资源量
    %求解这个线性规划问题,获得优化变量x
    %x中的每个元素表示对一个车辆任务使用对应的资源量
    %由于没有约束条件,x中的元素可取任意值,包含0和1之间的分数
    %然后后续代码根据x distinguish整数解和分数解,通过贪心方法获得近似最优解。
    %利用LP松弛是一种典型的求解组合优化问题的方法,可以快速提供一个可行解%
  x = linprog([], [], [], [], [], [], c, 0, [], optimoptions('linprog','Display', 'off'));
  
  J_1 = find(x==1);
  J_F = find(0<x<1);

  E_s_max = max(E_s(J_F)); 
  index = [];
  
  if sum(E_s(J_1)) >= E_s_max
    E_s = sum(E_s(J_1));
  else
    [E_s,index] = max(E_s(J_F));
  end
  
end
function h = channel_gain(d)

    alpha = 2; % 路径损耗指数
    d0 = 1; % 基准距离
    P_0 = 0.1; % 基准路径损耗
    
    h = P_0*(d0./d).^alpha;

end