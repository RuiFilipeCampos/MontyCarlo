import React, { useEffect, useState } from 'react'
import Link from 'next/link'
import Layout from '../components/Layout'
import { Text, Flex, VStack, HStack, Spacer } from '@chakra-ui/layout'
import { PlusSquareIcon, ViewIcon, ChevronDownIcon, ChevronRightIcon, MinusIcon, HamburgerIcon } from '@chakra-ui/icons'
import { Button, Icon, IconButton, Input, Textarea, useDisclosure } from '@chakra-ui/react'
import Editor  from "react-simple-code-editor";

import {
  Menu,
  MenuButton,
  MenuList,
  MenuItem,
  MenuItemOption,
  MenuGroup,
  MenuOptionGroup,
  MenuDivider,
} from '@chakra-ui/react'


import {
  Modal,
  ModalOverlay,
  ModalContent,
  ModalHeader,
  ModalFooter,
  ModalBody,
  ModalCloseButton,
} from '@chakra-ui/react'

function PopUp({children, buttonText}) {
  const { isOpen, onOpen, onClose } = useDisclosure()

  return (
    <>
      <Button onClick={onOpen}>{buttonText}</Button>

      <Modal closeOnOverlayClick={false} isOpen={isOpen} onClose={onClose}>
        <ModalOverlay />
        <ModalContent>
          <ModalHeader>Create a new folder</ModalHeader>
          <ModalCloseButton />
          <ModalBody pb={6}>
              {children}
          </ModalBody>

          <ModalFooter>
            <Button colorScheme='blue' mr={3}>
              Save
            </Button>
            <Button onClick={onClose}>Cancel</Button>
          </ModalFooter>
        </ModalContent>
      </Modal>
    </>
  )
}




const CreateFolderMenuItem = ({dirChildren, setDirChildren}) => {
  const { isOpen, onOpen, onClose } = useDisclosure()
  

  let [folderName, setFolderName] = useState('')
  const handleSubmit = () => {
      setDirChildren(
        [...dirChildren, <Directory name={folderName} />]
      )
      onClose()

  }

  return <>
  <MenuItem onClick={onOpen}>New Folder</MenuItem>
  <Modal closeOnOverlayClick={false} isOpen={isOpen} onClose={onClose}>
    <ModalOverlay />
    <ModalContent>
      <ModalHeader>Create a new folder</ModalHeader>
      <ModalCloseButton />
      <ModalBody pb={6}>
          <Input onChange={e=>setFolderName(e.target.value)}/>
      </ModalBody>

      <ModalFooter>
        <Button onClick={handleSubmit} colorScheme='blue' mr={3}>
          Save
        </Button>
        <Button onClick={onClose}>Cancel</Button>
      </ModalFooter>
    </ModalContent>
  </Modal>
</>
}

const FolderMenu = ({dirChildren, setDirChildren}) => {



  return <Menu>
    <MenuButton
        as={IconButton}
        aria-label='Options'
        icon={<HamburgerIcon />}
        variant='outline'
    />
    <MenuList>
      <MenuItem>New File</MenuItem>
      <CreateFolderMenuItem 
        dirChildren={dirChildren} 
        setDirChildren={setDirChildren}
      />
      <MenuItem>Delete</MenuItem>
    </MenuList>
  </Menu>

}




const File = ({selected=false}) => {
  if (selected) return <HStack align="center">
    <ViewIcon />
    <Text> filename.py </Text>
  </HStack>

  return <Button textAlign="left" variant="link">filename.py</Button>
}


const Directory = ({children, name, open=false}) => {
  let [isOpen, setIsOpen] = useState(open)

  const change = () => {
    if (isOpen) setIsOpen(false)
    else setIsOpen(true)
  }

  let [dirChildren, setDirChildren] = useState([<></>, <></>])


  if (isOpen) return <VStack align="left" w="full">
    <Flex align="center"  w="full">
      <VStack align="left" spacing={0.01}>
        <HStack align="left">
          <Button onClick={change} variant="link">{name}</Button>
          <Spacer/>
        </HStack>
        <HStack align="left">
          <Button size="xs" variant="link">New File</Button>
          <Button size="xs" variant="link">New Directory</Button>
        </HStack>
      </VStack>
    </Flex>
    <VStack align="left"  w="full">
      {children}
    </VStack>
  </VStack>


  return <VStack align="left" w="full">
    <Flex align="center" onClick={change} w="full">
    <VStack align="left" spacing={0.01}>
        <HStack align="left">
          <Button onClick={change} variant="link">{name}</Button>
          <Spacer/>
        </HStack>
        <HStack align="left">
          <Button size="xs" variant="link">New File</Button>
          <Button size="xs" variant="link">New Directory</Button>
        </HStack>
      </VStack>
    </Flex>
  </VStack>
}


const NewFolderButton = () => {
  let [folderName, setFolderName] = useState('')

  const onSave = () => {


  }

  return <PopUp buttonText="New Folder">
    <Input onChange={e=>setFolderName(e.target.value)}/> 
  </PopUp>
}




const DirectoryAndFiles = () => {

  let [jsx, setJsx] = useState([<></>, <></>])
  let [name, setName] = useState('')

  const createFolder = () => setJsx([...jsx, <Directory name={name}/>])

  return <div>
    <VStack align="left" h="full" w="full" >
      <Directory name="Your Project">
        <File />
      </Directory>
    </VStack>
  </div>
}


const OpenTerminal = () => <Flex bg="black" h="40%" w="full" overflow="hidden" >
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



const Terminal = ({terminalMode}) => {
if (terminalMode) return <OpenTerminal />

return <Flex bg="black" h="10%" w="full" overflow="hidden"  p={5} align="left" >
  <Text textColor="white">(your project name) C:\></Text>
</Flex>
}


function TextCode() {
  return <textarea style={{
      padding:15, 
      width:'99%', 
      height:"90%", 
      outline:'0px'
    }}>
      {'\n\n\ndef main(x):\n    return <geo.space fill="water">\n        <src.beam particle="photon" />\n    <geo.space/>'}
  </textarea>
}

const Coder = ({terminalMode}) => {
  if (terminalMode) return <Flex alignItems="center" h="50%" w="full" >
    <TextCode />
  </Flex>

  return <Flex alignItems="center" h="80%" w="full" bg="white">
    <TextCode />
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
  const toggleTerminalMode = () => {
    if (terminalMode) setTerminalMode(false)
    else setTerminalMode(true)
  }
  
  return <Layout title="Home | Next.js + TypeScript + Electron Example">
      <Flex h="inherit" w="20%">
        <VStack h="full" w="full" p={5} align="left" spacing={8}>
          <Text> The name of your project</Text>
          <DirectoryAndFiles />
        </VStack>
      </Flex>
      <VStack w="full" >
        <Coder terminalMode={terminalMode} />
        <Terminal terminalMode={terminalMode} />
        <HStack h="10%"  w="full" >
          <Button w="full" > Compile </Button>
          <Button w="full" > Plot Geometry </Button>
          <Button w="full" > Create new material </Button>
          <Button w="full" onClick={toggleTerminalMode} > Toggle Terminal </Button>
        </HStack>
      </VStack>
    </Layout>
}
export default IndexPage; 


