function fileNameList = getFileNameList(path)
    fileNameList = [];
    if isdir(path)
        fileSpecList = dir(path);
        fileNameListRaw = string({fileSpecList.name});
        fileNameListLen = length(fileNameListRaw);
        cleanListLen = fileNameListLen - ( any(strcmp(fileNameListRaw,".")) + any(strcmp(fileNameListRaw,"..")) );
        if cleanListLen==fileNameListLen
            fileNameList = fileNameListRaw;
        else
            fileNameList = strings([cleanListLen,1]);
            counter = 0;
            for i = 1:fileNameListLen
                if ~(strcmp(fileNameListRaw(i),".") || strcmp(fileNameListRaw(i),".."))
                    counter = counter + 1;
                    fileNameList(counter) = fileNameListRaw(i);
                end
            end
        end
    end
end