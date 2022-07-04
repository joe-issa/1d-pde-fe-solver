function Log (start_time)

txt_name = [start_time, '_log.txt'];

if ~exist("Logs", 'dir')
   mkdir("Logs")
end

cd Logs

diary(txt_name)

cd ..

end
