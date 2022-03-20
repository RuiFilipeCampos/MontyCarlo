import React, { useEffect, useState } from 'react'
import * as ch from '@chakra-ui/react'
import * as ch_icons from '@chakra-ui/icons'

import { FiSquare } from "react-icons/fi"
import {IoIosAddCircleOutline} from "react-icons/io"
import {AiFillFolderOpen} from 'react-icons/ai'
import {AiFillCaretRight, AiFillCaretDown} from 'react-icons/ai'











const originalFilesAndFolders = {
  name:'Top Level Folder',
  folders:[
    {
      name:'First Level Folder',
      folders:[],
      files:[ { pk:1, name:'filename.py', state:'closed', }, ], 
    },
  ],
  files:[ { pk:2, name:'filename.py', state:'closed', },  { pk:3, name:'filename.py', state:'closed', },],
}


const OpenedFiles = React.createContext(null)

import { createBreakpoints } from '@chakra-ui/theme-tools'
import { fileURLToPath } from 'url'



const print = (obj) => JSON.stringify(obj)

const breakpoints = createBreakpoints({
  sm: '320px',
  md: '768px',
  lg: '960px',
  xl: '1200px',
  '2xl': '1536px', 
})


var FILE_COUNTER = 0;

const getFileCounterValue = () => {
  FILE_COUNTER = FILE_COUNTER + 1
  return FILE_COUNTER
}

const File = ({this_file, its_folder}) => {
  const [forceUpdate, openedFiles, setOpenedFiles] = React.useContext(OpenedFiles)



  

  const minimize = () => {
    this_file.state = 'min'
  }

  this_file.minimize = minimize

  const handleClick = () => {
    switch (this_file.state) {

      case 'closed':
        for (let file of openedFiles) file.minimize()
        this_file.state = 'open' 
        setOpenedFiles( [ ...openedFiles, this_file])
        
        return

      
      case 'open':
        this_file.state = 'closed'

        let newOpenedFiles = []
        for (let file of openedFiles){
          if (file.pk == this_file.pk) continue
          newOpenedFiles.push(file)
        }

        if (newOpenedFiles.length > 0) newOpenedFiles[0].state = 'open'
        setOpenedFiles(newOpenedFiles)
        return

      case 'min':
        for (let file of openedFiles) file.minimize()
        
        this_file.state = 'open'
        forceUpdate()


        return
      

    
      default:
        break;
    }


  }



  if (this_file.state=='open') return <ch.Box paddingLeft={4} w="full">
    <ch.Text onClick={handleClick} roundedLeft="md" paddingLeft={4} bg='gray' textColor='black' px="1">
      {this_file.name}
    </ch.Text>
  </ch.Box>

  if (this_file.state=='min') return <ch.Box  paddingLeft={4} w="full">
    <ch.Text onClick={handleClick}  roundedLeft="md" paddingLeft={4} bg='whiteAlpha.300' textColor='black' px="1"         bg="blackAlpha.500" textColor="whiteAlpha.600"
>
      {this_file.name}
    </ch.Text>
  </ch.Box>

  return <ch.Box paddingLeft={4} w="full">
    <ch.Text onClick={handleClick} roundedLeft="md" textColor="whiteAlpha.500" _hover={{bg:'gray', textColor:'black', cursor:'pointer'}} px="1">
      {this_file.name}
    </ch.Text>
  </ch.Box>
}


const Folder = ({this_folder}) => {
  const [isOpen, setIsOpen] = React.useState(this_folder.isOpen)


  const handleClick = () => {
    this_folder.isOpen = !this_folder.isOpen

    setIsOpen(!isOpen)
  }


  if (!isOpen) return <ch.VStack paddingLeft={4}  align="left" w="full"> 
    <ch.HStack onClick={handleClick}   _hover={{bg:'whiteAlpha.200', textColor:'black', cursor:'pointer'}}>
      < AiFillCaretRight color="rgb(255, 255, 255, .3)"/>
      <ch.Text textColor="whiteAlpha.500">{this_folder.name}</ch.Text>
    </ch.HStack>
  </ch.VStack>



  let CONTENTS = []
  for (let folder of this_folder.folders){
    CONTENTS.push(
      <Folder this_folder={folder} />
    )
  }
  
  for (let file of this_folder.files){
    CONTENTS.push(
      <File 
        this_file={file}
        its_folder = {this_folder}
      />
    )
  }

  return <ch.VStack paddingLeft={4}  align="left" w="full"> 
      <ch.HStack onClick={handleClick}  _hover={{bg:'whiteAlpha.200', textColor:'black', cursor:'pointer'}}>
        < AiFillCaretDown color="rgb(255, 255, 255, .3)"/>
        <ch.Text textColor="whiteAlpha.500">{this_folder.name}</ch.Text>
      </ch.HStack>
      {CONTENTS}
  </ch.VStack>
}

 


const FileExplorer = ({filesAndFolders}) => {
  



  return <ch.VStack py={10} w="full" alignText="left" w="full">
      <Folder  this_folder = {filesAndFolders} />
  </ch.VStack> 
}


const TopBar = () =>     <ch.HStack overflow="clip" textOverflow="ellipsis" w="full" bg="black" px={2} style={{'-webkit-app-region': 'drag'}}>
<ch.Spacer />
<ch.Text overflow="clip" textOverflow="ellipsis"  textAlign="center" h="20px" fontSize="12px" w="full" textColor="whiteAlpha.500">
  Title of the Project
</ch.Text>
  <ch.Spacer />
  <ch.HStack spacing="1px">
    <ch.IconButton 
      style={{'-webkit-app-region': 'no-drag'}}
      _hover={{textColor:"white"}}
      variant="unstyled"
      textColor="whiteAlpha.500"
      icon={<ch_icons.MinusIcon/>}
      arial-label=""
    />
    <ch.Box 
      style={{'-webkit-app-region': 'no-drag'}}
      _hover={{textColor:"white"}}
      variant="unstyled"
      textColor="whiteAlpha.500"
      as={FiSquare}
      arial-label=""
    />
    <ch.IconButton 
      style={{'-webkit-app-region': 'no-drag'}}
      _hover={{'textColor':"red"}}
      variant="unstyled"
      textColor="whiteAlpha.500"
      icon={<ch_icons.CloseIcon/>}
      arial-label=""

    />
  </ch.HStack>
</ch.HStack>



const LeftBar = () => <ch.VStack 
  style={{'-webkit-app-region': 'drag'}}
  bg="black" h="full" w="50px" py={1}
  overflow="clip" spacing={5}
>
    <ch.Spacer/>
    <ch.Button 
      style={{'-webkit-app-region': 'no-drag'}}
      _hover={{textColor:"white", cursor:"pointer"}}
      variant="unstyled"
      textColor="whiteAlpha.500"
      as={ch_icons.PlusSquareIcon}
      arial-label=""
      size="xs"
    />
    <ch.Button 
      style={{'-webkit-app-region': 'no-drag'}}
      _hover={{textColor:"white", cursor:"pointer"}}
      variant="unstyled"
      textColor="whiteAlpha.500"
      as={AiFillFolderOpen}
      arial-label=""
      size="xs"
    />
    <ch.Button 
      style={{'-webkit-app-region': 'no-drag'}}
      _hover={{textColor:"white", cursor:"pointer"}}
      variant="unstyled"
      textColor="whiteAlpha.500"
      as={ch_icons.SearchIcon}
      arial-label=""
      size="xs"
    />
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.Spacer/>

    <ch.Spacer />
    <ch.IconButton 
        style={{'-webkit-app-region': 'no-drag'}}
        _hover={{textColor:"white"}}
        variant="unstyled"
        textColor="whiteAlpha.500"
        icon={<ch_icons.ChevronLeftIcon/>}
        arial-label=""
    />
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.Spacer/>
    <ch.IconButton 
      style={{'-webkit-app-region': 'no-drag'}}
      _hover={{textColor:"white"}}
      variant="unstyled"
      textColor="whiteAlpha.500"
      icon={<ch_icons.SettingsIcon/>}
      arial-label=""
    />
</ch.VStack>



const Line = ({index}) => {

  const [value, setValue] = React.useState('')
  let a = <ch.Flex h="1%" w="1%" bg="black"><>.</></ch.Flex>
  return  <ch.HStack bg="" w="full" h="25px" alignItems="center" spacing={4}>
    <ch.Text fontSize="sm" textColor="blackAlpha.400"> {index} </ch.Text>
    <ch.Box bg="" h="full" w='5000px'>
        <pre style={{
          'z-index':0,
          position: 'fixed',
          "font-size": "13pt",
          "font-family": "monospace",
          color:'black',
          overflow:'hidden',
          width:'5000px',
          height:"25px",
          'caret-color': 'black',
          cursor:'auto',
          background:'transparent',
        }}>{value}</pre>

        <ch.Input spellcheck="false"  w="5000px" outline="none" variant='unstyled' bg="transparent"
          style={{
            outline:"none", 
            border:'none',
            'z-index':1,
            overflow:'hidden',
            width:'5000px',
            // height:"100%",
            background:"transparent",
            position: 'absolute',
            color:'transparent',
            "font-size": "13pt",
            "font-family": "monospace",
            'caret-color': 'black',
            cursor:'pointer',
          }}
        
        onChange={(e)=>setValue(e.target.value)} 
        />



    </ch.Box>
  </ch.HStack>


}


const Body = () => {
  const [forceUpdate, openedFiles, setOpenedFiles] = React.useContext(OpenedFiles)
  const [value, setValue] = React.useState('')
  const handleTabsChange = index => {
    for (let file of openedFiles) file.minimize()
    openedFiles[index].state = 'open'
    forceUpdate()
  }

  let tabs = []
  let tabPanels = []
  let tabIndex = 0
  let i = 0

  for (let file of openedFiles) {


    if (file.state == 'open') tabIndex = i

    tabs.push(
      <ch.Tab 
        roundedTop="md" 
        _selected={{ textColor: 'black', bg: 'gray' }} 
        bg='whiteAlpha.300' textColor='black'
        bg="blackAlpha.500" textColor="whiteAlpha.600"

      >
        {file.name}
      </ch.Tab>
    )

    tabPanels.push(
      <ch.TabPanel w="full" h="full">
        <ch.VStack w="full" h="full" spacing="1px">
          <Line index={1} /> 
          <Line index={2} />
          <Line index={3} />
        </ch.VStack>
      </ch.TabPanel>
    )

    ++i
  }


  return <ch.VStack  h="full" w="full" overflow="hidden" textOverflow="clip"   >
    <ch.Tabs index={tabIndex} onChange={handleTabsChange} size='sm' variant='unstyled' w="full" bg="gray" h="full"   >
      <ch.TabList paddingTop={2}  bg="blackAlpha.600">
        {tabs}
      </ch.TabList>
      <ch.TabPanels h="full" w="full" >
        {tabPanels}
      </ch.TabPanels>
    </ch.Tabs>
    </ch.VStack>

}













const Index = () => {
  const [openedFiles, setOpenedFiles] = React.useState([])
  const [counter, setCounter] = React.useState(0)
  const forceUpdate = () => setCounter(counter + 1)


  return  <OpenedFiles.Provider value={ [forceUpdate, openedFiles, setOpenedFiles] }> 
    <ch.VStack w="100vw" h="100vh" spacing={0} >
      <TopBar />
      <ch.HStack bg="gray" w="full" h="full" spacing={0}>
        <LeftBar />
        <ch.VStack 
          bg="blackAlpha.600"
          h="full"
          w={["0px", "300px"]}
          overflow="hidden"
          textOverflow="clip"
        >
            <FileExplorer 
              filesAndFolders={originalFilesAndFolders}
            /> 
        </ch.VStack>
        <Body />
      </ch.HStack>
    </ch.VStack>
  </ OpenedFiles.Provider>

}

export default Index