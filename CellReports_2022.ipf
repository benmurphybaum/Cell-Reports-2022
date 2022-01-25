//Functions used for analysis in Murphy-Baum and Awatramani, 2022 in Cell Reports.

//Coded in IGOR PRO 8, although most of the functions are back-compatible with Igor 7

//These functions are meant to be used within the NeuroTools interface (DOI in the paper) and the functions won't compile without it.
//However, the important mathematical and logical operations used should be apparent and reproducible.



#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function NT_PeakTime(DS_ROIs,PercentPeak,FitTime,StartTime,EndTime)
	String DS_ROIs //Data set containing Ca signal waves.
	Variable PercentPeak //0 to 1, 1 just takes the location peak itself
	Variable FitTime //amount of time prior to the peak location to begin the sigmoid fit
	Variable StartTime,EndTime
	
	//SUBMENU=Turning Project
	//TITLE=Peak Time
	
//	Note={
//	Finds the time point of the peak value, or some fraction of the peak value.
//	Performs a sigmoid fit to the peak event to more accurately get the peak rise times.
//	Output is a single wave with the peak times for each wave set.
//	
//	\f01PercentPeak\f00 : Time point when the amplitude is this fraction of the actual peak. 0-1.
//	      -Set to 1 to take the actual peak time. 
//	\f01FitTime\f00 : Amount of time prior to the peak to start the sigmoid fit (1-2 s)
//	\f01StartTime\f00 : Start X value for finding the peak
//	\f01EndTime\f00 : Ending X value for finding the peak
//	}
//	
	STRUCT ds ds
	GetStruct(ds)
	
	ds.wsi = 0
	
	If(PercentPeak > 1 || PercentPeak < 0)
		Abort "PercentPeak must be between 0 and 1"
	EndIf
	
	SetDataFolder GetWavesDataFolder(ds.waves[0],1)
	
	String outName = "PkLoc_" + NameOfWave(ds.waves[0])
	Make/O/N=(ds.numWaves[0]) $outName/Wave=outWave
	
	
	Note outWave,"Percent Peak: " + num2str(PercentPeak)
	Note outWave,"Peak Location Waves:"
	
	Do
		Wave ROI = ds.waves[ds.wsi]
		SetDataFolder GetWavesDataFolder(ROI,1)
		
		WaveStats/Q/R=(StartTime,EndTime) ROI
		
		//peak location X scale
		Variable pk = V_MaxLoc
		
		//Peak amplitude
		Variable pkVal = V_Max
		
		//peak location X point
		Variable pkPnt = ScaleToIndex(ROI,pk,0)
		
		Variable fitStart = pk - FitTime //fit starts 1 second prior to the peak 
		Variable fitStartPt = ScaleToIndex(ROI,fitStart,0)
		
		//Fit a sigmoid to the data for reducing measurement noise
		try
			CurveFit/Q sigmoid ROI[fitStartPt,pkPnt]/D;AbortOnRTE
		catch
			Variable error = GetRTError(1)
			print GetErrMessage(error),"...continuing..."
			continue
		endtry	
		
		Wave fit = $("fit_" + NameOfWave(ROI))
		
		If(!WaveExists(fit))
			Abort "Couldn't find the fit wave"
		EndIf
		
		Variable threshold = PercentPeak * pkVal
		FindLevel/EDGE=1/Q fit,threshold
		
		If(V_flag)
			print "Couldn't find the level. Using NaN"
			outWave[ds.wsi] = nan
		Else
			outWave[ds.wsi] = V_LevelX
		EndIf
		
		KillWaves/Z fit
		
		Note outWave,NameOfWave(ROI)
			
		ds.wsi += 1
	While(ds.wsi < ds.numWaves[0])
	
End

//Counts spikes for turning stimuli - Change Direction Project
Function NT_TurnSpikeCount(DS_SpikeData,Threshold,StartTime,EndTime,FolderName,menu_Range,menu_TurnType,filterTerm)
	//SUBMENU=Turning Project
	//TITLE=Spike Count (Turn Data)
	
//	Note={
//	Counts spikes, names the output waves specifically for turning stimuli.
//	
//	\f01Range\f00 : Specifies time range for the spike count
//	\f01TurnType\f00 : Specifies linear, ± 90°, or 180° turns
//	\f01FilterTerm\f00 : String for additional filtering of a data set. Waves that
//	    don't match the filter will be ignored.
//	}
	
	//Data Set for the spiking data, raw spiking data waves
	String DS_SpikeData
	
	//Some input variables
	Variable Threshold,StartTime,EndTime
	
	//Name of the folder to put the spike counts in
	String FolderName
	
	//Drop down menu for the time range we're counting in - for wave naming purposes
	String menu_Range
	
	//Drop down menu for the type of turn the data contains
	String menu_TurnType
	
	//This will filter out certain wavesets. If the string match fails on the first wave, it will skip that wave set
	String filterTerm
	
	String menu_Range_List = "early;late;all;"
	String menu_TurnType_List = "Linear;Turn90;TurnN90;Turn180;"
	
	String Threshold_Assign = "root:Packages:NeuroToolsPlus:ControlWaves:threshold"
	
	//Data set info structure
	STRUCT ds ds 
	
	//Fills the data set structure
	GetStruct(ds)
	
	//Reset wave set index
	ds.wsi = 0
	
	//Filter term match test
	If(strlen(filterTerm))
		If(!stringmatch(NameOfWave(ds.waves[0]),filterTerm))
			return 0
		EndIf
	EndIf
	
	String folder = GetWavesDataFolder(ds.waves[0],1)
	
	If(strlen(FolderName))
		folder = ParseFilePath(1,folder,":",1,0) // back out one folder
	
		folder += FolderName
	EndIf
	
	If(!DataFolderExists(folder))
		NewDataFolder $folder
	EndIf
	
	SetDataFolder $folder
	
	//Name of the output wave that will hold the results
	String outputName = "DSSpk_" + menu_Range + "_" + StringsFromList("2-*",NameOfWave(ds.waves[0]),"_",noEnding=1)
	Make/O/N=(ds.numWaves[0]) $outputName/Wave = outWave
	
	String vAngName = "vAng_" + menu_Range + "_" + StringsFromList("2-*",NameOfWave(ds.waves[0]),"_",noEnding=1)
	Make/O/N=(1) $vAngName/Wave = vAng
	
	String dsiName = "vDSI_" + menu_Range + "_" + StringsFromList("2-*",NameOfWave(ds.waves[0]),"_",noEnding=1)
	Make/O/N=(1) $dsiName/Wave = dsi
	
	//Function Loop
	Do
		If(endTime == 0)
			endTime = pnt2x(ds.waves[ds.wsi],DimSize(ds.waves[ds.wsi],0)-1)
		EndIf
		
		//Spike count
		outWave[ds.wsi] = GetSpikeCount(ds.waves[ds.wsi],StartTime,EndTime,Threshold)
		
		ds.wsi += 1
	While(ds.wsi < ds.numWaves[0])
	
	//Were any spikes detected? If not, toss out the tuning curve
	If(sum(outWave) == 0)
		KillWaves/Z vAng,dsi,outWave
		return 0
	EndIf
	
	SetScale/P x,0,45,"deg",outWave
	
	//Vector Sum angle
	vAng[0] = VectorSum(outWave,"0;45;90;135;180;225;270;315;","angle")
	//Vector Sum DSI
	dsi[0] = VectorSum(outWave,"0;45;90;135;180;225;270;315;","DSI")
	
	
	Note outWave,"Threshold: " + num2str(Threshold)
	Note outWave,"StartTime: " + num2str(StartTime)
	Note outWave,"EndTime: " + num2str(EndTime)
	Note outWave,"Range: " + menu_Range
	Note outWave,"TurnType: " + menu_TurnType
	Note outWave,"PD: " + num2str(vAng[0])
	Note outWave,"DSI: " + num2str(dsi[0])
End

//Counts spikes for turning stimuli - Change Direction Project
Function NT_CollisionSpikeCount(DS_SpikeData,Threshold,StartTime,EndTime,menu_Range)
	//SUBMENU=Turning Project
	//TITLE=Spike Count (Collision Data)
	
//	Note={
//	Counts spikes, names the output waves specifically for colliding stimuli.
//	
//	\f01Range\f00 : Specifies time range for the spike count
	
	//Data Set for the spiking data, raw spiking data waves
	String DS_SpikeData
	
	//Some input variables
	Variable Threshold,StartTime,EndTime

	//Drop down menu for the time range we're counting in - for wave naming purposes
	String menu_Range
	
	String menu_Range_List = "early;collision;late;"
	
	String Threshold_Assign = "root:Packages:NeuroToolsPlus:ControlWaves:threshold"
	
	//Data set info structure
	STRUCT ds ds 
	
	//Fills the data set structure
	GetStruct(ds)
	
	//Reset wave set index
	ds.wsi = 0
	
	String folder = GetWavesDataFolder(ds.waves[0],1)
	SetDataFolder $folder
	
	//Name of the output wave that will hold the results
	String outputName = "DSSpk_" + StringFromList(1,NameOfWave(ds.waves[0]),"_") + "_" + menu_Range + "_" + StringsFromList("2-*",NameOfWave(ds.waves[0]),"_",noEnding=1)
	Make/O/N=(ds.numWaves[0]) $outputName/Wave = outWave
	
	AddOutput(outWave,ds)
	
	//Function Loop
	Do
		If(endTime == 0)
			endTime = pnt2x(ds.waves[ds.wsi],DimSize(ds.waves[ds.wsi],0)-1)
		EndIf
		
		//Spike count
		outWave[ds.wsi] = GetSpikeCount(ds.waves[ds.wsi],StartTime,EndTime,Threshold)
		
		ds.wsi += 1
	While(ds.wsi < ds.numWaves[0])
		
	Note outWave,"Threshold: " + num2str(Threshold)
	Note outWave,"StartTime: " + num2str(StartTime)
	Note outWave,"EndTime: " + num2str(EndTime)
	Note outWave,"Range: " + menu_Range
End


//Calculates the difference in angle between two waves. 
Function NT_AngleDistance(DS_Turns,menu_Range,Suffix)
	//SUBMENU=Turning Project
	//TITLE=Angular Distance
	
//	Note = {
//	Gets the angular distance between the apparent preferred directions for -90° 
//	and +90° turns, compared with the preferred direction for the linear stimulus.
//	
//	Each wave set must contain all three (+90°,-90°, and linear) angle varieties.
//	i.e. the function takes three waves at a time.
//
//	\f01Range\f00 : Time range that the tuning curves were acquired in (naming purposes)
//	\f01Suffix\f00 : Suffix applied to the output wave, defaults to 'delta'
//	}
	
	String DS_Turns //Data set containing linear, turn90, and turnN90 avg angles in each wave set.
	//If all three aren't present, we ignore the wave set.
	
	String menu_Range,Suffix
	String menu_Range_List = "early;late;all;"
	
	//Data set info structure
	STRUCT ds ds 
	
	//Fills the data set structure
	GetStruct(ds)
	
	//Data set checks, must have three waves
	If(ds.numWaves[0] != 3)
		return 0
	EndIf
	
	Variable i
	
	If(!strlen(suffix))
		suffix = "delta"
	EndIf
	
	//allocate to waves
	For(i=0;i<3;i+=1)
		String name = NameOfWave(ds.waves[i])
		
		If(stringmatch(name,"*t90*"))
			Wave turn90 = ds.waves[i]
			String outName = name + "_" + suffix
			SetDataFolder GetWavesDataFolder(turn90,1)
			Duplicate/O turn90,$outName
			Wave turn90_delta = $outName
			
		ElseIf(stringmatch(name,"*t270*"))
			Wave turnN90 = ds.waves[i]
			outName = name + "_" + suffix
			SetDataFolder GetWavesDataFolder(turnN90,1)
			Duplicate/O turnN90,$outName
			Wave turnN90_delta = $outName
			
		Else
			Wave linear = ds.waves[i]
		EndIf
	EndFor
	
	//Compute the signed difference between them
	turn90_delta[0] = polarMath(linear[0],turn90[0],"deg","distance",1)
	
	turnN90_delta[0] = polarMath(linear[0],turnN90[0],"deg","distance",1)
	
End

Function NT_CircularStats(DS_Waves,cb_Radians,cb_outputInDegrees)
	//SUBMENU=Turning Project
	//TITLE=Circular Stats
	
	//Returns the circular mean, median, and stdev of a distribution of angles
	
	String DS_Waves //input waves containing angular data
	Variable cb_Radians //is the data already in radians? 1 or 0 (i.e. True or False)
	Variable cb_outputInDegrees //do you want the output data to be in degrees? 1 or 0 (i.e. True or False)
	
	//Data set info structure
	STRUCT ds ds 
	
	//Fills the data set structure
	GetStruct(ds)
	
	Variable i
	For(i=0;i<ds.numWaves[0];i+=1)
		Wave theWave = ds.waves[i]
		
		SetDataFolder GetWavesDataFolder(theWave,1)
		
		//Make radians scaled wave for the mean angles
		If(!cb_Radians)
			Duplicate/O theWave,$(ds.paths[i][0] + "_rad")
			Wave rad = $(ds.paths[i][0] + "_rad")
			rad = theWave * pi/180
		Else
			Wave rad = theWave
		EndIf		
		
		//Extra column needed for vector lengths, set all to 1
		Redimension/N=(-1,2) rad
		rad[][1] = 1
		
		StatsCircularMeans/Z/CI rad
		Wave stats = W_CircularMeans
		
		Duplicate/O stats,$(ds.paths[i][0] + "_CircularMeans")
		KillWaves/Z stats
		
		Wave stats = $(ds.paths[i][0] + "_CircularMeans")
		If(cb_outputInDegrees)
			stats[1] = stats[1] * 180/pi
			stats[2] = stats[2] * 180/pi
			stats[3] = stats[3] * 180/pi
			stats[4] = stats[4] * 180/pi
		EndIf
		
		Redimension/N=(-1,1) rad
		rad = (rad < 0) ? rad + 2*pi : rad
		StatsCircularMoments/Q/Z rad
		Wave stats = W_CircularStats
		
		If(cb_outputInDegrees)
		
			stats[10] = stats[10] * 180/pi//circular standard deviation
			stats[8] = stats[8] * 180/pi //mean
			stats[11] = stats[11] * 180/pi //median
		EndIf
		
		Duplicate/O stats,$(ds.paths[i][0] + "_CircularStats")
		KillWaves/Z stats

	EndFor
End

//Calculates the modulation index of two Ca signals. Calculates the peak first, then a % difference.
//Used in Figure 3 of the manuscript
Function NT_Modulation_Index2(DS_Data1,DS_Data2,StartTime,EndTime,PeakWidth,cb_Abs_Value,Output_Suffix)
	//SUBMENU=Turning Project
	//TITLE=Modulation Index 2
//	
//	Note={
//	Calculates the modulation index between the peaks of two waves (difference / sum)
//	Corresponding waves in each data set are compared. Same as Modulation Index, but uses
// two data sets instead of one to define the data. 
//	
//	\f01StartTime:\f00 Starting X value for finding the peak
//	\f01EndTime:\f00 Ending X value for finding the peak
//	\f01PeakWidth:\f00 Size of window to average around the peak value.
//	}
	
	//DS_Data1 is a data set holding responses to Trial 1 of the turning stimulus
	//DS_Data2 is a data set holding responses to Trial 2 of the turning stimulus (reversed path)
	
	String DS_Data1,DS_Data2
	Variable StartTime,EndTime,PeakWidth,cb_Abs_Value
	String Output_Suffix
	
	//Data set info structure
	STRUCT ds ds 
	
	//Fills the data set structure
	GetStruct(ds)
	
	//Reset wave set index
	ds.wsi = 0
	
	DFREF saveDF = GetDataFolderDFR()
	
	If(ds.numWaves[0] != ds.numWaves[1])
		Abort "Each waveset must have the same number of waves"
	EndIf
	
	//modulation index wave is put in the first wave's folder
	Wave wave1 = ds.waves[0][0] //data 1
	SetDataFolder GetWavesDataFolder(wave1,1)
	
	If(strlen(Output_Suffix))
		Output_Suffix = "_" + Output_Suffix
	EndIf
	
	Make/O/N=(ds.numWaves[0]) $("MI" + Output_Suffix)/Wave=MI
	Note/K MI,"Modulation Index (difference / sum)"
	Note MI,"Waves:"
	
	Variable pk1,pk2,modulationIndex
	
	PeakWidth /= 2
		
	Do
		//declare each wave in the wave set
		Wave wave1 = ds.waves[ds.wsi][0] //data 1
		Wave wave2 = ds.waves[ds.wsi][1] //data 2
		
		//Data 1
		WaveStats/Q/R=(StartTime,EndTime) wave1 //find the peak location
		
		//Get the average value ±PeakWidth around the peak location
		pk1 = mean(wave1,V_maxLoc - PeakWidth,V_maxLoc + PeakWidth)
		
		//Data 2
		WaveStats/Q/R=(StartTime,EndTime) wave2 //find the peak location
		
		//Get the average value ±PeakWidth around the peak location
		pk2 = mean(wave2,V_maxLoc - PeakWidth,V_maxLoc + PeakWidth)
		
		//Ensure no negative peak values
		pk1 = (pk1 < 0) ? 0 : pk1
		pk2 = (pk2 < 0) ? 0 : pk2
		
		//calculate modulation index for absolute value (unsigned DSI) or signed DSI.
		If(cb_Abs_Value)
			MI[ds.wsi] = abs(pk1 - pk2) / (pk1 + pk2)
		Else
			MI[ds.wsi] = (pk1 - pk2) / (pk1 + pk2)
		EndIf

		Note MI,NameOfWave(wave1) + " vs. " + NameOfWave(wave2)
		
		ds.wsi += 1
	While(ds.wsi < ds.numWaves[0])
	
	SetDataFolder saveDF
	
End


//Used in Figure 3 of the manuscript to determine ROI position relative to the location of the stimulus turn
Function NT_DefineROISector(DS_BaseImage,CDF_MI,X_Center,Y_Center,menu_Split,menu_ROI_Group,menu_diagonalAngle)
	//SUBMENU=Turning Project
	//TITLE=Define ROI Sector
	
//	Note = {
//	For a given set of ROIs, it assigns each ROI to different spatial sectors,
//	as defined by the 'Split' input, and the X/Y centers. Modulation index
//	for these ROIs is sorted into one wave for each spatial sector.
//	
//	\f01MI\f00 : Modulation index wave. Must be same size as number of ROIs
//	\f01X_Center\f00 : X center point on the image from where the divide is drawn
//	\f01Y_Center\f00 : Y center point on the image from where the divide is drawn
//	\f01Split\f00 : Defines how the image is split up into sectors
//	\f01ROI_Group\f00 : ROI Group from previously defined ROIs
//	\f01diagonalAngle\f00 : Defines angle to divide the image if diagonal split is chosen
//	}
	
	String DS_BaseImage,CDF_MI
	Variable X_Center,Y_Center
	String menu_Split //Should we split the receptive field vertically, horizontally, or diagonally
	String menu_ROI_Group
	String menu_diagonalAngle
	
	String menu_Split_List = "Vertical;Horizontal;Diagonal;"
	String menu_diagonalAngle_List = "45;-45;" //+45 is bottom left to top right, -45 is top left to bottom right
	
	//Use the ROI lists as the menu items 
	String menu_ROI_Group_List = TextWaveToStringList(root:Packages:NeuroToolsPlus:ScanImage:ROIGroupListWave,";")
	
	Variable angle = str2num(menu_diagonalAngle)
	
	DFREF NTSI = $SI
	
	//Data set info structure
	STRUCT ds ds 
	
	//Fills the data set structure
	GetStruct(ds)

	//If the data set happens to have more than one image in it, just use the first one
	Wave theWave = ds.waves[0]
	
	If(!WaveExists(theWave))
		return 0
	EndIf
	
	//Make sure it's a 2D wave
	If(WaveDims(theWave) < 2)
		return 0
	EndIf
	
	//Make sure the MI wave is selected and valid
	Wave/Z MI = $CDF_MI
	If(!WaveExists(MI))
		print "Couldn't find the MI Wave"
		return 0
	EndIf	
	
	//YOUR CODE GOES HERE....
	
	Variable rows =  DimSize(theWave,0)
	Variable cols =  DimSize(theWave,1)
	
	//Resolve the X center point
	If(X_Center > 1e-3 && X_Center < 1) //above scale of the image (mm not microns) means its a fractional input
		Variable xTurn = IndexToScale(theWave,DimSize(theWave,0) * X_Center,0)
	Else
		//Make sure center point is actually within the range of the images
		Variable index = ScaleToIndex(theWave,X_Center,0)
		If(index > rows || index < 0)
			DoAlert 1,"Turn center point was not within the image scale. Continue anyway?"
			
			If(V_flag == 2)
				Abort ""
			EndIf
		EndIf
		
		xTurn = X_Center
	EndIf
	
	//Resolve the Y center point
	If(Y_Center > 1e-3 && Y_Center < 1) //above scale of the image (mm not microns) means its a fractional input
		Variable yTurn = IndexToScale(theWave,DimSize(theWave,1) * Y_Center,1)
	Else
		//Make sure center point is actually within the range of the images
		index = ScaleToIndex(theWave,Y_Center,1)
		If(index > cols || index < 0)
			DoAlert 1,"Turn center point was not within the image scale. Continue anyway?"
			
			If(V_flag == 2)
				Abort ""
			EndIf
		EndIf
	
		yTurn = Y_Center
	EndIf
	
	//Get the center XY coordinates of the ROI group
	SI_GetCenter(group =menu_ROI_Group)
	
	SVAR software = NTSI:imagingSoftware
	
	strswitch(software)
		case "2PLSM":
			Wave xROI = $("root:twoP_ROIs:" + menu_ROI_Group + "_ROIx")
			Wave yROI = $("root:twoP_ROIs:" + menu_ROI_Group + "_ROIy")
			break
		case "ScanImage":
			Wave xROI = $("root:Packages:NeuroToolsPlus:ScanImage:ROIs:" + menu_ROI_Group + "_ROIx")
			Wave yROI = $("root:Packages:NeuroToolsPlus:ScanImage:ROIs:" + menu_ROI_Group + "_ROIy")
			break
	endswitch
	
	
	If(!WaveExists(xROI) || !WaveExists(yROI))
		return 0
	EndIf
	
	//Assign each ROI to one side of the image or the other
	Variable numROIs = DimSize(xROI,0)
	If(numROIs == 0)
		return 0
	EndIf
	
	//Name of the output wave that will hold the results
	SetDataFolder root:Analysis:$menu_ROI_Group
	String outputName = menu_ROI_Group + "_Sectors"
	
	DFREF MIFolder = GetWavesDataFolderDFR(MI)
	
	//Make the output wave 
	Make/O/N=(numROIs) $outputName/Wave = outWave
	
	Variable i,count1=0,count2=0
	For(i=0;i<numROIs;i+=1)
		Variable ycoord = yROI[i]
		Variable xcoord = xROI[i]
		
		strswitch(menu_Split)
			case "Vertical":
				If(i == 0)
					Make/O/N=1 :top/Wave=out1
					Make/O/N=1 :bottom/Wave=out2
				EndIf
						
				If(ycoord > yTurn)
					outWave[i] = 1 //top
					out1 += MI[i]
					count1 += 1
				Else
					outWave[i] = 0 //bottom
					out2 += MI[i]
					count2 += 1
				EndIf
				
				//Split the original MI wave into sectors for scatter plotting
				Duplicate/O MI,MIFolder:MI_top
				Duplicate/O MI,MIFolder:MI_bottom
				
				Wave MI_1 = MIFolder:MI_top
				Wave MI_2 = MIFolder:MI_bottom
				break
			case "Horizontal":
				If(i == 0)
					Make/O/N=1 :left/Wave=out1
					Make/O/N=1 :right/Wave=out2
				EndIf
				
				If(xcoord < xTurn)
					outWave[i] = 1 //left
					out1 += MI[i]
					count1 += 1
				Else
					outWave[i] = 0 //right
					out2 += MI[i]
					count2 += 1
				EndIf
				
				//Split the original MI wave into sectors for scatter plotting
				Duplicate/O MI,MIFolder:MI_left
				Duplicate/O MI,MIFolder:MI_right
				
				Wave MI_1 = MIFolder:MI_left
				Wave MI_2 = MIFolder:MI_right
				
				break
			case "Diagonal":
				If(i == 0)
					Make/O/N=1 :above/Wave=out1
					Make/O/N=1 :below/Wave=out2
				EndIf
				
				Variable slope = sin(angle * pi/180) / cos(angle * pi/180)
				Make/FREE line
				SetScale/I x,-200e-6,200e-6,line //set the scale range ±200 microns
				
				line =  (x - X_Center) * slope + Y_Center //y = mx + b
				
				//check if ROI is below the line
				Variable threshold = line[x2pnt(line,xCoord)]
				
				If(yCoord > threshold)
					//above the line
					outWave[i] = 1
					out1 += MI[i]
					count1 += 1
					
				Else
					//below the line
					outWave[i] = 0
					out2 += MI[i]
					count2 += 1
				EndIf
				
				//Split the original MI wave into sectors for scatter plotting
				Duplicate/O MI,MIFolder:MI_above
				Duplicate/O MI,MIFolder:MI_below
				
				Wave MI_1 = MIFolder:MI_above
				Wave MI_2 = MIFolder:MI_below
				
				break
		endswitch
	EndFor

	out1 /= count1
	out2 /= count2
	
	Variable num1 = sum(outWave) //number of sites in left/top/above
	Variable num2 = DimSize(outWave,0) - num1 //number of sites in right/bottom/below
	Redimension/N=(num1) MI_1
	Redimension/N=(num2) MI_2
	
	count1 = 0
	count2 = 0
	
	MI_1 = 0
	MI_2 = 0
	
	For(i=0;i<DimSize(outWave,0);i+=1)
		If(outWave[i] == 1)
			MI_1[count1] = MI[i] //left or top or above
			count1 += 1
		Else
			MI_2[count2] = MI[i] //right or bottom or below
			count2 += 1
		EndIf
	EndFor
	
	KillWaves/Z left,right
End


Function NT_DistanceFromSoma(DS_ROI_X,DS_ROI_Y,Soma_X,Soma_Y)
	String DS_ROI_X,DS_ROI_Y //X and Y coordinate waves for all the ROIs. Should be single waves.
	Variable Soma_X,Soma_Y //Scaled coordinates of the soma or whatever reference point we're using
	
	//SUBMENU=Turning Project
	//TITLE=Distance From Soma
	
	//	Note={
	//	Calculates the distance (line of sight, not cable distance) of ROIs from a target location,
	//	usually the soma but could be any coordinate. 
	//	
	//	Input two waves, one with all the ROIs X coordinates and the other with the Y coordinates.
	//	}
		
	STRUCT ds ds
	GetStruct(ds)
	
	If(ds.numWaves[%ROI_X] != 1 || ds.numWaves[%ROI_Y] != 1 )
		DoAlert 0,"Must have single X ROI coordinates wave and a single Y ROI coordinates wave."
		return 0
	EndIf
	
	//Get the coordinates waves
	Wave xROI = ds.waves[0][%ROI_X]
	Wave yROI = ds.waves[0][%ROI_Y]
	
	//Make the output distance wave
	String folder = GetWavesDataFolder(xROI,1)
	
	String outPath =  folder + RemoveEnding(NameOfWave(xROI),"x") + "_DistanceToSoma"
	Make/O/N=(DimSize(xROI,0)) $outPath/Wave=distance
	
	//calculate distance from the target coordinate.
	distance = sqrt((xROI[p] - Soma_X)^2 + (yROI[p] - Soma_Y)^2)

End

//Counts the number of spikes in the wave
//This is for other functions to use, operates on a single wave
Function GetSpikeCount(theWave,StartTime,EndTime,Threshold)
	Wave theWave
	Variable StartTime,EndTime,Threshold
	
	If(!WaveExists(theWave))
		return -1
	EndIf
	
	DFREF saveDF = GetDataFolderDFR()
	SetDataFolder GetWavesDataFolder(theWave,1)
	
	Duplicate/FREE theWave,temp
	
	FlattenWave(temp)
	
	FindLevels/Q/D=spktm/R=(StartTime,EndTime)/M=0.002/T=0.0005 temp,threshold

	Variable spkct = V_LevelsFound
	KillWaves/Z W_FindLevels,spktm
	
	SetDataFolder saveDF
	
	return spkct
End

//Returns the vector sum angle or dsi of the input wave
Function [Variable PD,Variable DSI,Variable Radius] VectorSum2(Wave theWave,String angles)
	
	Variable i,size = DimSize(theWave,0)
	
	If(ItemsInList(angles,";") != size)
		return [-1,-1,-1] //angle list must be the same length as the input wave
	EndIf
	
	Variable vSumX,vSumY,totalSignal
	
	vSumX = 0
	vSumY = 0
	totalSignal = 0

	For(i=0;i<size;i+=1)
		If(numtype(theWave[i]) == 2) //protects against NaNs, returns -9999, invalid
			return [-9999,-9999,-9999]
		EndIf
		
		Variable theAngle = str2num(StringFromList(i,angles,";"))
		
		vSumX += theWave[i]*cos(theAngle*pi/180)
		vSumY += theWave[i]*sin(theAngle*pi/180)
		totalSignal += theWave[i]
	EndFor
	
	Variable vRadius = sqrt(vSumX^2 + vSumY^2)
	Variable vAngle = atan2(vSumY,vSumX)*180/pi
	DSI = vRadius/totalSignal
	Variable SNR = vRadius
	
	If(vAngle < 0)
		vAngle +=360
	Endif

	Make/N=3/O $(GetWavesDataFolder(theWave,2) + "_VectorSum")  /Wave=VSData
	VSData[0] = vAngle
	VSData[1] = DSI
	VSData[2] = vRadius
	
	SetDimLabel 0,0,Angle,VSData
	SetDimLabel 0,1,DSI,VSData
	SetDimLabel 0,2,Radius,VSData
	
	PD = vAngle
	Radius = vRadius
	
	return [PD,DSI,Radius]
End

//Get ROI ------------------------------------
//Extracts the time-varying fluorescence from selected scans at specified ROIs.
Function/WAVE NT_GetROI(menu_Measure,menu_Channel,BaselineStart,BaselineEnd,menu_FilterType,Filter,bt_Display_Last)
	String menu_Measure,menu_Channel //Measurement type (∆F/F in the paper), channel #
	Variable BaselineStart,BaselineEnd //Start and End times for the baseline region
	String menu_FilterType //type of filter to temporally smooth the Ca data, S-G 2nd in the paper (Savitsky Golay 2nd order)
	Variable Filter //Filter magnitude, 9 in the paper
	String bt_Display_Last //Auto display the extracted data
	
	String menu_Measure_List = "∆F/F;SDEV;ABS;"
	String menu_Channel_List = "1;2;1/2;2/1;"
	String menu_FilterType_List = "None;S-G 2nd;S-G 4th;Boxcar;Gaussian;"
	
	String bt_Display_Last_Pos = "30;220;120;20"
	String bt_Display_Last_Proc = "DisplayLastROIsProc"
	
	//SUBMENU=Imaging
	//TITLE=Get ROI
	
//	Note={
//	Extracts the time-varying fluorescence from selected scans at specified ROIs.
//	
//	Use the Image Browser to select scans and the ROIs before running (Ctrl-1)
//	
//	\f01Measure\f00 : Measurement type (∆F/F, Standard Dev. above baseline, or absolute F)
//	\f01Channel\f00 : Channel of the scan, and defines ∆G/R etc. potentially
//	\f01BaselineStart\f00 : Start of baseline region for ∆F calculation, in seconds
//	\f01BaselineEnd\f00 : End of baseline region for ∆F calculation, in seconds
//	\f01Filter\f00 : Savitsky-Golay filter, odd numbered greater than 5
//	     -9 is a good starting point usually.
//	}
	
	STRUCT ds ds
	GetStruct(ds)
	
	STRUCT IMAGING img
	
	initParam(img,ds)
	
	DFREF NPC = $CW
	DFREF NTSI = $SI
	DFREF NTD = $DSF
	
	DFREF saveDF = GetDataFolderDFR()

	//Make the ROI analysis folder if it doesn't already exist
	If(!DataFolderExists("root:Analysis"))
		NewDataFolder root:Analysis
	EndIf
	SetDataFolder root:Analysis	
	
	//Make the output wave reference wave for passing the result onto another function
	Make/FREE/WAVE/N=0 outputWaveRefs
	Make/O/WAVE/N=(img.scan.num * img.roi.num) NTSI:roiRefs/Wave=roiRefs
	
	Variable i,j,k,totalWaveCount = 0

	//SCAN LOOP
	For(i=0;i<img.scan.num;i+=1)
		Variable ref = StartMSTimer
		
		strswitch(menu_Channel)
			case "1": //channel 1 only
				Wave theScan = img.scan.ch1[i] //signal fluorescence
				Wave theBgnd = img.scan.ch1[i] //background fluorescence
				break
			case "2": //channel 2 only
				Wave theScan = img.scan.ch2[i]
				Wave theBgnd = img.scan.ch2[i]
				break
			case "1/2": // ch1 / ch2
				Wave theScan = img.scan.ch1[i]
				Wave theBgnd = img.scan.ch2[i]
				break
			case "2/1": // ch2 / ch1
				Wave theScan = img.scan.ch2[i]
				Wave theBgnd = img.scan.ch1[i] 
				break
		endswitch

		//Get dendritic mask
		Wave mask = GetDendriticmask(theBgnd)
		Redimension/B/U mask
		
		//Get dark value
		ImageStats/R=mask theBgnd
		Variable darkVal = 0.9*V_avg
		
		//Cleanup
		KillWaves/Z mask,root:Analysis:maxProj
		
		//ROI LOOP
		For(j=0;j<img.roi.num;j+=1)
			String theROI = img.rois[j][0][0]
			String roiGroup = img.rois[j][0][1]
			
			//Make the ROI Group Analysis folder
			If(!DataFolderExists("root:Analysis:" + roiGroup))
				NewDataFolder $("root:Analysis:" + roiGroup)
			EndIf
			
			//Make the ROI analysis subfolder
			String ROIFolder = "root:Analysis:" + roiGroup + ":" + theROI
			
			If(!DataFolderExists(ROIFolder))
				NewDataFolder $ROIFolder
			EndIf
			
			//X and Y waves that define the ROI area
			Wave roiX = img.roi.x[j]
			Wave roiY  = img.roi.y[j]
			
			If(DimSize(roiX,1) > 0)
				Variable isSoma = 1
			Else
				isSoma = 0
			EndIf
			
			//Use somatic ROI map instead of the traditional boundary ROIs
			If(isSoma)
				SetDataFolder $ROIFolder
//				NT_GetROI_Soma(roiX,theScan,theBgnd,img)
				SetDataFolder saveDF
				continue
			EndIf
			
			//Seed values for filling out the ROI mask
			Variable maskMax,maskMin,xSeed,ySeed
			WaveStats/Q theBgnd
			
			maskMin = WaveMin(roiX)
			maskMax = WaveMax(roiX)
			
			xSeed = maskMax + DimDelta(theBgnd,0)
			If(xSeed > IndexToScale(theBgnd,DimSize(theBgnd,0)-1,0))
				xSeed = IndexToScale(theBgnd,0,0)
			EndIf
			
			maskMin = WaveMin(roiY)
			maskMax = WaveMax(roiY)
			
			ySeed = maskMax + DimDelta(theBgnd,1)
			If(ySeed > IndexToScale(theBgnd,DimSize(theBgnd,1)-1,1))
				ySeed = IndexToScale(theBgnd,0,1)
			EndIf
			
			//ROI mask wave	
			SetDataFolder $ROIFolder			
			ImageBoundaryToMask ywave=roiY,xwave=roiX,width=(DimSize(theBgnd,0)),height=(DimSize(theBgnd,1)),scalingwave=theBgnd,seedx=xSeed,seedy=ySeed			
		
			Wave ROIMask = $(ROIFolder + ":M_ROIMask")	
			
			//Did the ROI mask actually get created?
			If(!WaveExists(ROIMask))
				DoAlert 0, "Couldn't find the ROI mask wave for: " + NameOfWave(theScan)
				continue
			EndIf
			
			//Make the raw ROI waves for signal and background
			Variable numFrames = DimSize(theScan,2)
			Make/O/FREE/N=(numFrames) ROI_Signal,ROI_Bgnd
			
			//Set all the scales of the ROI waves
			SetScale/P x,DimOffset(theScan,2),DimDelta(theScan,2),ROI_Signal,ROI_Bgnd
			
			//Average values over the ROI region
			For(k=0;k<numFrames;k+=1)
				ImageStats/M=1/P=(k)/R=ROImask theScan
				ROI_Signal[k] = V_avg
				
				ImageStats/M=1/P=(k)/R=ROImask theBgnd
				ROI_Bgnd[k] = V_avg
			EndFor		
			
			//Temporal smoothing options
			If(Filter)
				strswitch(menu_FilterType)
					case "None":
						break
					case "S-G 2nd":
						If(Filter < 5)
							print "2nd Order Savitsky-Golay filter requires odd values >= 5. Setting to 5."
							Filter = 5
						EndIf
						Smooth/S=2 (Filter), ROI_Signal,ROI_Bgnd
						break
					case "S-G 4th":
						If(Filter < 7)
							print "4th Order Savitsky-Golay filter requires odd values >= 7. Setting to 7."
							Filter = 7
						EndIf
						Smooth/S=4 (Filter), ROI_Signal,ROI_Bgnd
						break
					case "Boxcar":
						Smooth/B (Filter), ROI_Signal,ROI_Bgnd
						break
					case "Gaussian":
						Smooth (Filter), ROI_Signal,ROI_Bgnd
						break
				endswitch
			EndIf
					
			//Use median for the baseline, so it doesn't get pulled up or down from noisy values
			Variable	bsln = median(ROI_Bgnd,BaselineStart,BaselineEnd)
			
			//Absolute fluorescence or delta fluorescence?
			strswitch(menu_Measure)
				case "∆F/F":
					//∆F/F
					String outName = NameOfWave(theScan) + "_" + theROI + "_dF"
					break
				case "SDEV":
					//Standard Deviation
			 		outName = NameOfWave(theScan) + "_" + theROI + "_sd"
					break
				case "ABS":
					//Abs
					outName = NameOfWave(theScan) + "_" + theROI + "_abs"
					break
			endswitch
		
			//Make the dF or dG wave
			Make/O/N=(numFrames) $outName
			Wave dF = $outName
			
			//Set all the scales of the ROI waves
			SetScale/P x,DimOffset(theScan,2),DimDelta(theScan,2),dF
			
			//Calculate the ∆F/F or Absolute fluoresence ratios
			
			darkVal = 0
			
			strswitch(menu_Measure)
				case "∆F/F":
					//∆F/F
					dF = (ROI_Signal - bsln) / (bsln - darkVal)
					break
				case "SDEV":
					//Standard Deviation
			 		
			 		//baseline subtracted and dark subtracted signal
					dF = ROI_Signal - bsln - darkVal
					
					//get standard deviation of the baseline region
					WaveStats/Q/R=(BaselineStart,BaselineEnd) dF 
					
					//standard deviations above the median baseline value
					dF /= V_sdev
					break
				case "ABS":
					//Abs
					dF = ROI_Signal
					break
			endswitch
			
			//Set the wave note with the original scan name
			Note dF,"Scan: " + GetWavesDataFolder(theScan,2)
			Note dF,"TYPE: " + menu_Measure
			Note dF,"FILTER TYPE: " + menu_FilterType
			Note dF,"FILTER SIZE: " + num2str(Filter)
			Note dF,"BASELINE START: " + num2str(BaselineStart)
			Note dF,"BASELINE END: " + num2str(BaselineEnd)
			
			//These are all the output ROI waves
			Redimension/N=(totalWaveCount + 1) outputWaveRefs
			outputWaveRefs[totalWaveCount] = dF
			roiRefs[totalWaveCount] = dF
			
			totalWaveCount += 1
		EndFor

		print "Get ROI (" +  roiGroup + "):",NameOfWave(theScan) + ",",StopMSTimer(ref) / (1e6),"s"
	EndFor
	
	SetDataFolder saveDF
		
	//pass the output wave on
	return outputWaveRefs
End

//inputs spike recording, outputs histograms
Function/WAVE NT_PSTH(DS_Waves,BinSize,SpikeThreshold,menu_Type,OutputFolder,cb_RemoveTrend,StartTime,EndTime)
	//SUBMENU=Spiking
	
//	Note={
//	Calculates the spike rate over time. Input is raw spiking data.
//	
//	\f01BinSize:\f00 Size of the bin to measure the firing rate, 0.02 (20 ms) is good default.
//	\f01SpikeThreshold:\f00 Threshold value for detecting spikes. Make sure units are correct.
//	\f01Type:\f00 Either binned or convolved with a gaussian kernel of width = BinSize
//	\f01OutpuFolder:\f00 Output folder for the PSTH, put inside the folder with the data
//	\f01RemoveTrend:\f00 Eliminates wobbly baseline using a polynomial fit
//	\f01StartTime:\f00 Starting X point for the PSTH
//	\f01EndTime:\f00 Ending X point for the PSTH. 0 will take the entire wave
//	}
	
	String DS_Waves //raw spiking data
	Variable BinSize,SpikeThreshold
	String menu_Type,OutputFolder //menu type is Gaussian or Binned PSTH calculation, Gaussian used in the paper
	Variable cb_RemoveTrend,StartTime,EndTime
	
	String menu_Type_List = "Gaussian;Binned;"
	String SpikeThreshold_Assign = "root:Packages:NeuroToolsPlus:ControlWaves:threshold"
	String StartTime_Assign = "root:Packages:NeuroToolsPlus:ControlWaves:rangeLeft"
	String EndTime_Assign = "root:Packages:NeuroToolsPlus:ControlWaves:rangeRight"
	
	STRUCT ds ds
	GetStruct(ds)
	
	Variable i,j,numWaves,numBins
	
	SetDataFolder GetWavesDataFolder(ds.waves[0],1)
	
	//Check start and end time validity
	If(EndTime == 0 || EndTime < StartTime)
		EndTime = pnt2x(ds.waves[0],DimSize(ds.waves[0],0) -1)
	EndIf
	
	
	//Reset the wsi
	ds.wsi = 0
	Do
		Wave theWave = ds.waves[ds.wsi]
		
		//Which data folder should we put the output waves?
		If(strlen(OutputFolder))
			CreateFolder(GetWavesDataFolder(ds.waves[ds.wsi],1) + OutputFolder)
		EndIf
			
		SetDataFolder GetWavesDataFolder(ds.waves[ds.wsi],1) + OutputFolder
		
		//Remove low pass trends in the wave to flatten it
		If(cb_RemoveTrend)
			FlattenWave(theWave)
		EndIf
		
		//Make Spike Count wave
		Make/FREE/N=(ds.numWaves[0]) spkct
		
		//Get spike times and counts
		FindLevels/Q/EDGE=1/M=0.002/R=(StartTime,EndTime)/D=spktm theWave,SpikeThreshold
		spkct[i] = V_LevelsFound
		
		//Gaussian or binned histograms
		strswitch(menu_Type)
			case "Binned":	
				numBins = floor((IndexToScale(theWave,DimSize(theWave,0)-1,0) - IndexToScale(theWave,0,0) )/ BinSize) //number of bins in wave
				String histName = RemoveEnding(ReplaceListItem(0,NameOfWave(theWave),"_","PSTH"),"_")
				Make/O/N=(numBins) $histName
				Wave hist = $histName
				
				//add to output data set
				AddOutput(hist,ds)
				
				If(DimSize(spktm,0) == 0)
					hist = 0
				Else
					Histogram/C/B={pnt2x(theWave,0),BinSize,numBins} spktm,hist
				EndIf
				
				hist /= binSize
				
				break
			case "Gaussian":
				Variable dT = DimDelta(theWave,0)
				Variable sampleRate = 1000 // 1 ms time resolution
				//gaussian template for convolution
				Make/O/N=(3*(BinSize*sampleRate)+1) template
				Wave template = template
				SetScale/I x,-1.5*BinSize,1.5*BinSize,template
				template = exp((-x^2/(0.5*BinSize)^2)/2)
				
				Variable theSum = sum(template)
				template /= (1000*binSize)
				
				Variable histDelta = (DimSize(theWave,0)*dT)/sampleRate
				Make/O/FREE/N=(DimSize(theWave,0)*dT*sampleRate) raster
				
				SetScale/P x,0,1/sampleRate,raster
				raster = 0
				
				For(j=0;j<DimSize(spktm,0);j+=1)
					If(x2pnt(raster,spktm[j]) > (DimSize(raster,0)-1))
						continue
					EndIf
					raster[x2pnt(raster,spktm[j])] = 1
				Endfor
				
				histName = RemoveEnding(ReplaceListItem(0,NameOfWave(theWave),"_","PSTH"),"_")
				Duplicate/O template,$histName
				Wave hist = $histName
				
				//add to output data set
				AddOutput(hist,ds)
				
				Convolve raster, hist
				hist *=1000
				
				break
		endswitch	
		
		//Cleanup
		KillWaves spktm,template
		
		//Set the wave note
		String noteStr = "PSTH:\r"
		noteStr += "Type: " + menu_Type + "\r"
		noteStr += "Threshold: " + num2str(SpikeThreshold) + "\r"
		noteStr += "Bin Size: " + num2str(BinSize) + "\r"
		noteStr += "StartTm: " + num2str(StartTime) + "\r"
		noteStr += "EndTm: " + num2str(EndTime) + "\r"
		
		Note/K hist,noteStr
				
		ds.wsi += 1
	While(ds.wsi < ds.numWaves[0])
	
End