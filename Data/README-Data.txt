In the Data folder, a CSV file can be found with the hydrochemical data obtained and analysed in this article. 

The header columns of the CSV file describe the data which can be observed below. Extra information about the type of data in these columns can be found below:

Column A: Sample Text - Sample identity as given during the fieldwork (type: string)
Column B: #Number - Sample number (type: int)
Column C: Name - This column shows in which monitoring well the sample is taken. The names used are slightly different than in the article. 
		PZ1 = MW3
		PZ2 = MW6
		PZ3 = MW1
		PZ4 = MW4
		PZ5 = MW5
		PZ6 = MW2
Column D: Date - Sampling date (type: date dd-mm-yyyy)
Column E: Time - Sampling time (type: time hh:mm)
Column F: Code - Codename given to a specific type of sample, often GW (groundwater) or INF (infiltrated water), but during the PPT more different codes are used like: DW=Tile Drainage Water, DW+=Tile Drainage Water including additions, SP=Sample (type: string)
Column G: Type - Type of sample, can be Sample or Sensor. Sample is data obtained from a watersample and sensor is obtained sensor data (type: string)
Column H: Project - Data was obtained from different projects like PPT (Push-Pull Test), BST1 (CFIE), BST2 (storage period - Winter 2019), BST4 (storage period - Fall 2020), BST5 (storage period - Spring 2021).
Column I: O2 (mg/l) - O2 concentrations obtained from Ponsel OPTOD sensor (type: float)
Column J: NO3 (mg/l) - NO3 concentrations obtained from TRiOS OPUS or NICO NO3 sensor. NOTE: data not always reliable! (type: float)
Column K: DOC (mg/l) - DOC concentrations obtained from TRiOS OPUS. NOTE: data not always reliable! (type: float)
Column L: TOC (mg/l) - DOC concentrations obtained from TOC-V CPH analyser, Shimadzu, Japan  (type: float)
Column M: T (C) - Temperature obtained from sensor in degrees Celcius. (type: float)
Column N: pH - pH values obtained from Ponsel PHEHT sensor (type: float)
Column O: eH (mV) - eH values obtained from Ponsel PHEHT sensor (type: float)
Column P: EC (uS/cm) - Electrical conductivity in uS/cm obtained from Ponsel C4E sensor (type: float)
Column Q: Alk (pulses/ 5ml) - Manual titration for alkanity concentration, number of pulses needed before color change in 5 ml watersample (only used during PPT) (type: int)
Column R: Abs. Volume - Only used during PPT, abstracted volume during this push-pull test time step in L
Column S: Total abs. Volume - Only used during PPT, total abstracted volume during this push-pull test time step plus the previous in L
Column T-AW: Elemental analysis in watersample by Inductively Coupled Plasma – Mass Spectrometry (ICP-MS; PlasmaQuant MS, Analytik-Jena, Germany), concentrations in ug/L (type: float)
Column AX-BB: Analysis of cations by Ion Chromatography, concentrations in mg/L (type: float)
Column BC-BJ: Analysis of anions by Ion Chromatography, concentrations in mg/L (type: float)
Column BK: Elemental analyis of P by Inductively Coupled Plasma - optical emission spectrometry, concentrations in mg/L (type: float)
Column BL-CK: Analysis of pesticides by Liquid Chromatography – Mass Spectrometry (LC-MS; Xevo TQ-S micro, Waters, U.S.A.), concentrations in ug/L (type: float)
Column CL-CN: Analysis of NH4, PO4, Alkalinity obtained from Discrete Analysis (DA; AQ400, Seal analytical, UK), concentrations in mg/L (type: float)  
