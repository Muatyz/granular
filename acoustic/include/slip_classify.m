function type=slip_classify(shear_stress_drop, slip_event_label, method, threshould)
%%% To classify the slip events into micro-, minor- and major-slips.
%%%
%%% `shear_stress_drop`: the record of the ss drop (in time sequence);
%%% `slip_event_label`: used to classify the ss drop data, such as cut-off
%%% frequency, frequency peak, 
%%% `method`:   - "ssDrop"
%%%             - "peakFreq"
%%%             - "cutoffFreq"
%%% `threshould`: micro- | minor- | major- 
%%% `type`: 0 | 1 | 2

num = length(shear_stress_drop);
type = zeros(1, num);

if method == "ssDrop"
    %%% 此时 threshould 应为观察应力降直方图/时间分布图后人为区分的
    %%% 应力区间, 如 threshould = [0.2, 0.4] 分为 micro\in (0,0.2],
    %%% minor\in (0.2, 0.4], major(failure)\in (0.4,\infty)
    %%% 在该分类下, `shear_stress_drop`=`slip_event_label`; 其余情况下
    %%% 则是维度相同.
    for i = 1:num
        if slip_event_label(i) <= threshould(1)
            type(i) = 0;
        elseif slip_event_label(i) > threshould(1) && slip_event_label(i) <= threshould(2)
            type(i) = 1;
        else
            type(i) = 2;
        end
    end
else
    error("Classifier function is not completed. Try other parametes to finish the process.")
end

end