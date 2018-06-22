%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rsPCA toolbox for helium-pump artifact removal
%
% Deaprtment of Brain and Cognitive Engineering, Korea University 
% Brain Signal Processing Laboraty,BSPL
%
% updated 07/22/2014
%
% Any suggestions or errors, please contact us, hyunchul_kim@korea.ac.kr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Please make sure that EEGLAB has been installed and EEG data is loaded before
% conducting rsPCA
%
% Example:
%     >> addpath ..\rsPCA -end
%     >> addpath ..\EEGLAB -end
%     >> rspca(EEG) % load an EEG data via EEGLAB
%
%

function rspca(EEG)

fs = round(EEG.srate);
sigpdetc = [0.01 0.02 0.03];
result_f = [];
flg_verbose = 0; % verbose off

f = figure('MenuBar', 'None','Visible','off','position',[360, 500, 450,285]);

segsize =[fs*2 fs*3 fs*4];

stgsigpdetc = ['1%'; '2%'; '3%'];
stgsegsize = [sprintf('%04d / 2s',fs*2); sprintf('%04d / 3s',fs*3); sprintf('%04d / 4s',fs*4)];
%% default parameters

TR = 2;  % default TR
seg_val = round(EEG.srate*2); % EEG segment of samples
sigp_val = 0.01;     % percentage level for the single peak detection
output_path = EEG.filepath; % output directory

%%

htext_seg  = uicontrol('Style','text','String','EEG segment size/second(s)',...
    'Position',[240,160,90,40]);

hpopup_segsize = uicontrol('Style','popupmenu','String',stgsegsize, ...
    'Position',[245,135,80,20],'Callback',{@popup_segsize_Callback});

htext_sigp  = uicontrol('Style','text','String','Percentage threhold level(%)',...
    'Position',[340,160,100,40]);

hpopup_sigp_detc = uicontrol('Style','popupmenu','String',stgsigpdetc, ...
    'Position',[355,130,70,25],'Callback',{@popup_sigp_Callback});

hrun = uicontrol('Style','pushbutton','String','Run', ...
    'Position',[340,50,100,50],'Callback',{@popup_run_Callback});

hresult = uicontrol('Style','pushbutton','String','Result', ...
    'Position',[340,10,100,30],'Callback',{@popup_result_Callback});

htext_list_box = uicontrol('style','text','String','EEG electrodes', ...
    'Position',[10,220,100,20]);

htext_box_input  = uicontrol('style','text','String','rsPCA inupt', ...
    'Position',[130,220,100,20]);

hlist_box = uicontrol('style','listbox','String',cellstr({EEG.chanlocs.labels}), ...
    'Position',[10,10,100,200],'Callback',{@hlist_box_Callback});

hlist_box_input = uicontrol('style','listbox','String',[], ...
    'Position',[130,10,100,200],'Callback',{@hlist_box_input_Callback});

hall_box = uicontrol('Style','pushbutton','String','ALL', ...
    'Position',[260,40,50,20],'Callback',{@popup_all_box_Callback});

hreset_box = uicontrol('Style','pushbutton','String','Reset', ...
    'Position',[260,20,50,20],'Callback',{@popup_all_reset_box_Callback});

hcheck_verbose_box = uicontrol('Style','checkbox','String',[], ...
    'Position',[280,75,15,15],'Callback',{@popup_check_verbose_box_Callback});

htext_verbose_checkbox  = uicontrol('style','text','String','Verbose', ...
    'Position',[250,95,70,15]);

htext_fs = uicontrol('style','text','String',sprintf('fs: %2.2fHz',fs), ...
    'Position',[220,220,120,20]);

htext_tr_txt = uicontrol('style','text','String','TR:', ...
    'Position',[340,220,30,20]);

htext_tr =uicontrol('style','edit','String',num2str(TR), ...
    'Position',[370,220,60,20]);

hset_output_path_box = uicontrol('Style','pushbutton','String','Output directory', ...
    'Position',[10,250,140,30],'Callback',{@popup_output_path_box_Callback}); 

htext_output_path_box = uicontrol('Style','text','String',output_path, ...
    'Position',[160,250,280,30]);

% Initialize the GUI.
% Change units to normalized so components resize automatically.

group_box = [hpopup_segsize,hpopup_sigp_detc, htext_sigp, htext_list_box, htext_box_input...
    htext_output_path_box, hset_output_path_box, hrun,htext_seg,hlist_box ...
    , hresult, hlist_box_input,hall_box,hreset_box, hcheck_verbose_box, ...
    htext_fs, htext_verbose_checkbox,htext_tr_txt,htext_tr];


set([f,group_box],'Units','normalized');

% Assign the GUI a name to appear in the window title.
set(f,'Numbertitle','off','Name','rsPCA for helium-pump artifact removal')

% Move the GUI to the center of the screen.
movegui(f,'center');

% Make the GUI visible.
set(f,'Visible','on');

%  Pop-up menu callback.
    function popup_segsize_Callback(source,eventdata)
        % Determine the selected data set.
        val = get(source,'Value');
        % Set current data to the selected data set.
        seg_val = segsize(val);
        disp(seg_val)
        disp(sprintf('EEG segement size = %d', segsize(val)));
        %
    end

    function popup_sigp_Callback(source,eventdata)
        val2 = get(source,'Value');
        sigp_val = sigpdetc(val2);
        disp(sprintf('the pecetrage threshold level = %2.2d %%', sigpdetc(val2)*100));
%         disp('call_sigp');
    end

    function hlist_box_Callback(source,eventdata)
        list_entry = cellstr(get(source,'String'));
        index_selected = get(source,'Value');  %changed line
        choice_listbox1 = list_entry(index_selected);
        update_listbox2 = cellstr(get(hlist_box_input, 'String'));
        
        flg = strcmp(choice_listbox1,update_listbox2);
        
        if length(update_listbox2)==1 && isempty(cell2mat(update_listbox2))
            newmenu = choice_listbox1;
            set(hlist_box_input,'String', newmenu);
        elseif flg==0
            newmenu = [update_listbox2 ;choice_listbox1 ];
            set(hlist_box_input,'String', newmenu);
        end
        
    end

    function hlist_box_input_Callback(source,eventdata)
        currentItems = cellstr(get(source,'String'));
        index_selected = get(source,'Value');  %changed line
        newItems =    currentItems;
        
        newItems(index_selected) = [];
        if isempty(newItems) == 1
            set(hlist_box_input,'String','');
        elseif strcmp(newItems,'');
            set(source,'Value',1)
            set(hlist_box_input,'String','');
        else
            set(source,'Value',length(newItems))
            set(hlist_box_input,'String',newItems);
        end
    end

    function popup_output_path_box_Callback(source,eventdata)
%         currentItems = cellstr(get(source,'String'));
        output_path = uigetdir;
        htext_output_path_box = uicontrol('Style','text','String',output_path, ...
            'Position',[160,250,280,30]);
        set([f,hpopup_segsize,hpopup_sigp_detc, htext_sigp, htext_list_box, htext_box_input...
            hset_output_path_box,htext_output_path_box, hrun,htext_seg,hlist_box, hlist_box_input,hall_box, hreset_box],'Units','normalized');
    end

    function popup_all_box_Callback(source,eventdata)
        list_entry = cellstr(get(hlist_box,'String'));
        set(hlist_box_input,'String',list_entry);
    end

    function popup_all_reset_box_Callback(source,eventdata)
        list_entry = cellstr(get(hlist_box,'String'));
        set(hlist_box_input,'Value',1);
        set(hlist_box_input,'String','');
    end

    function popup_check_verbose_box_Callback(source,eventdata)
        flg_verbose = get(source,'value');
        %         disp(list_entry);
    end

    function popup_run_Callback(source,eventdata)
        list_entry = cellstr(get(hlist_box_input,'String'));
        
        if strcmp(list_entry,'')
            h=errordlg('There is no input channel!','Error');
            set(h, 'WindowStyle', 'modal');
            uiwait(h);
        else
            
            %% find target channel 
            tgch = [];
            for list_idx=1:length(list_entry)
                for ch_idx=1:length(EEG.chanlocs)
                    if strcmp(list_entry(list_idx),EEG.chanlocs(ch_idx).labels);
                        tgch = [tgch ch_idx];
                    end
                end
            end
            
            stgch = sort(tgch);
            
            set(group_box,'enable','off');
            disp('run rsPCA!!');
            %             h = waitbar(0,'Please wait...');

            main_rspca(EEG,outdir,tgch,1,seg_val,sigp_val,true,flg_verbose)

            set(group_box,'enable','on');
            
        end
    end

    function popup_result_Callback(source,eventdata)
        TR=str2num(get(htext_tr,'string'));
        rspca_result_gui(EEG,TR,output_path);
    end

end


