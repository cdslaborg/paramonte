function arch = getArch()
    computerArch = computer("arch");
    if contains(computerArch,"64")
        arch = "x64";
    elseif contains(computerArch,"32")
        arch = "x32";
    else
        arch = computerArch;
    end
end