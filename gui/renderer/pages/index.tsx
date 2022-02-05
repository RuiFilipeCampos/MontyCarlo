import React, { useEffect, useState } from 'react'
import Link from 'next/link'
import Layout from '../components/Layout'
import { Text, Flex, VStack, HStack, Spacer } from '@chakra-ui/layout'
import { PlusSquareIcon, ViewIcon, ChevronDownIcon, ChevronRightIcon } from '@chakra-ui/icons'

const File = ({selected=false}) => {
  if (selected) return <HStack align="center">
    <ViewIcon />
    <Text> filename.py </Text>
  </HStack>

  return <Text _hover={{bg:'gray'}}> filename.py</Text>
}


const Directory = ({children, open=false}) => {
  let [isOpen, setIsOpen] = useState(open)

  const change = () => {
    if (isOpen) setIsOpen(false)
    else setIsOpen(true)
  }

  if (isOpen) return <VStack align="left" w="full">
    <Flex align="center" onClick={change} w="full">
      <ChevronDownIcon />
      <Text >Name of Directory</Text>
      <Spacer/>
    </Flex>
    <VStack align="left" px={4} w="full">
      {children}
    </VStack>
  </VStack>


  return <VStack align="left" w="full">
    <Flex align="center" onClick={change} w="full">
      <ChevronRightIcon />
      <Text >Name of Directory</Text>
      <Spacer/>
    </Flex>
  </VStack>
}

const DirectoryAndFiles = () => <div style={{ overflow: "scroll", height: "full" }}>
<VStack align="left" h="full" w="full" >
  <Directory open={true}>
    <File selected={true}/>
    <File />
    <File />
  </Directory>
  <Directory> 
    <File />
  </Directory>
  <Directory open={true}>
    <Directory open={true}> 
      <File />
    </Directory>
    <File />
  </Directory>
  <File />
  <File />
</VStack>
</div>



const OpenTerminal = () => <Flex bg="black" h="800vh" w="full" overflow="hidden" >
  <VStack p={5} align="left">
    <Spacer/>
    <Text textColor="white">(your project name) C:\> ls</Text>
    <Text textColor="white">
    {"CHANGELOG.md STYLEGUIDE.md config do.py manage.py reports"}
    </Text>
    <Text textColor="white">(your project name) C:\></Text>
    <Text textColor="white">(your project name) C:\></Text>
    <Text textColor="white">(your project name) C:\></Text>
  </VStack>
</Flex>


const Terminal = ({terminalMode, onClick}) => {
if (terminalMode) return <OpenTerminal />

return <Flex bg="black" h="100vh" w="full" onClick={onClick} overflow="hidden" >
  <Text textColor="white">(env) C:\></Text>
</Flex>
}

const Coder = ({onClick}) => {
  return <Flex onClick={onClick} h="full" w="full">
    <Text>Code</Text>
  </Flex>
}


const IndexPage = () => {
  useEffect(() => {
    // add a listener to 'message' channel
    global.ipcRenderer.addListener('message', (_event, args) => {
      alert(args)
    })
  }, [])



  let [terminalMode, setTerminalMode] = useState(false)



  return <Layout title="Home | Next.js + TypeScript + Electron Example">
      <Flex h="inherit" w="20%" shadow="lg">
        <VStack h="full" w="full" p={5} align="left" spacing={8}>
          <Text> The name of your project</Text>
          <DirectoryAndFiles />
        </VStack>
      </Flex>
      <VStack w="full"  >
        <Coder onClick={()=>setTerminalMode(false)}/>
        <Spacer/>
        <Terminal terminalMode={terminalMode} onClick={()=>setTerminalMode(true)}/>
      </VStack>
    </Layout>
}

export default IndexPage
