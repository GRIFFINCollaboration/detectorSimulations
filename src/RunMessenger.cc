///////////////////////////////////////////////////////////////////////////
//// The Messenger
//#include "RunMessenger.hh"
//#include "RunAction.hh"
//#include "G4UIdirectory.hh"
//#include "G4UIcmdWithABool.hh"
//#include "G4UIcmdWithAnInteger.hh"
//#include "G4UIcmdWithAString.hh"

//RunActionMessenger::RunActionMessenger(RunAction* pTarget)
//:myTarget(pTarget)
//{ 
//  myDirectory = new G4UIdirectory("/DetSys/run/");
//  myDirectory->SetGuidance("Control of run parameters.");

//  RunNumberCmd = new G4UIcmdWithAnInteger("/DetSys/run/runNumber", this);
//  RunNumberCmd->SetGuidance("Set the run number.");
//  RunNumberCmd->SetGuidance("Required parameters: 1 integer.");
//  RunNumberCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
//    
//  EnableWriteCmd = new G4UIcmdWithABool("/DetSys/run/enableLM", this);
//  EnableWriteCmd->SetGuidance("Activate the list-mode output.");
//  EnableWriteCmd->SetGuidance("Required parameters: none.");
//  EnableWriteCmd->SetParameterName("writeLMD",true);
//  EnableWriteCmd->SetDefaultValue(true);
//  EnableWriteCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

//  DisableWriteCmd = new G4UIcmdWithABool("/DetSys/run/disableLM", this);
//  DisableWriteCmd->SetGuidance("Deactivate the list-mode output.");
//  DisableWriteCmd->SetGuidance("Required parameters: none.");
//  DisableWriteCmd->SetParameterName("writeLMD",true);
//  DisableWriteCmd->SetDefaultValue(false);
//  DisableWriteCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

//  EnableWriteHexCmd = new G4UIcmdWithABool("/DetSys/run/enableHEX", this);
//  EnableWriteHexCmd->SetGuidance("Activate the list-mode output.");
//  EnableWriteHexCmd->SetGuidance("Required parameters: none.");
//  EnableWriteHexCmd->SetParameterName("writeHEX",true);
//  EnableWriteHexCmd->SetDefaultValue(true);
//  EnableWriteHexCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

//  DisableWriteHexCmd = new G4UIcmdWithABool("/DetSys/run/disableHEX", this);
//  DisableWriteHexCmd->SetGuidance("Deactivate the list-mode output.");
//  DisableWriteHexCmd->SetGuidance("Required parameters: none.");
//  DisableWriteHexCmd->SetParameterName("writeHEX",true);
//  DisableWriteHexCmd->SetDefaultValue(false);
//  DisableWriteHexCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

//  outputTimeUnit = new G4UIcmdWithAString("/DetSys/run/timeUnit",this);
//  outputTimeUnit->SetGuidance("Set the time output units, ie. second, millisecond, etc.");
//  outputTimeUnit->SetParameterName("choice",false);
//  outputTimeUnit->AvailableForStates(G4State_PreInit,G4State_Idle);
//}

//RunActionMessenger::~RunActionMessenger()
//{
//  delete myDirectory;
//  delete EnableWriteCmd;
//  delete DisableWriteCmd;
//  delete EnableWriteHexCmd;
//  delete DisableWriteHexCmd;
//  delete RunNumberCmd;
//  delete outputTimeUnit;
//}

//void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
//{ 
//  if( command == EnableWriteCmd ) {
//    myTarget->EnableWrite( EnableWriteCmd->GetNewBoolValue(newValue) );
//  }
//  if( command == DisableWriteCmd ) {
//    myTarget->EnableWrite( DisableWriteCmd->GetNewBoolValue(newValue) );
//  }
//  if( command == EnableWriteHexCmd ) {
//    myTarget->EnableWriteHex( EnableWriteHexCmd->GetNewBoolValue(newValue) );
//  }
//  if( command == DisableWriteHexCmd ) {
//    myTarget->EnableWriteHex( DisableWriteHexCmd->GetNewBoolValue(newValue) );
//  }
//  if( command == RunNumberCmd ) {
//    myTarget->SetRunNumber( RunNumberCmd->GetNewIntValue(newValue) );
//  }
//  if( command == outputTimeUnit ) {
//    myTarget->SetTimeUnit(newValue);
//  }
//}


